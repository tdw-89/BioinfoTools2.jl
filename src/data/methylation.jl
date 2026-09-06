"""
Single-base methylation calls.

A per-sample signal like the rest of [`Data`](@ref), but a purely *positional*
one: nothing here touches a `Genome`, a `Scaffold` or any of `Data`'s interval
types, and the only thing shared with the rest of the package is the `BitCodes`
packing vocabulary. Keep it that way, so a whole-genome call set stays one flat,
8-byte-per-site, memory-mappable block.
"""
module Methylation

using Arrow
using CodecZlib
using StructArrays

# Nested one level deeper than the other modules, so `BitCodes` is three dots up.
using ...BitCodes

#= Bit layout =#

# Field layout of the 32-bit call payload (see `AggregatedCall`). Offsets are
# 0-based bit positions.
const DEPTH_SHIFT = 0
const DEPTH_WIDTH = 16
const PERCENT_SHIFT = 16
const PERCENT_WIDTH = 8
const CONTEXT_SHIFT = 24
const CONTEXT_WIDTH = 2
const STRAND_SHIFT = 26
const STRAND_FIELD_WIDTH = STRAND_WIDTH   # 2 bits; bits 28-31 stay reserved

"""
Largest total read depth the 16-bit depth field can hold. Depths above this
saturate rather than wrapping into a neighbouring field; the methylation level
is encoded before the clamp and so stays accurate either way.
"""
const MAX_COUNT = UInt32(65535)

"""
Largest value of the 8-bit methylation-percentage field, i.e. the code for 100%
methylated. The field is a plain linear map of `[0, 100]` percent onto
`[0, 255]` — see [`encode_percent`](@ref) / [`decode_percent`](@ref).
"""
const MAX_PERCENT_CODE = UInt8(255)

"""
Cytosine sequence contexts, as stored in the 2-bit `context` field of an
[`AggregatedCall`](@ref) payload.

| Code   | Context                         |
|--------|---------------------------------|
| `0x00` | CpG                             |
| `0x01` | CHG                             |
| `0x02` | CHH                             |
| `0x03` | unknown (Bismark `U`/`u` calls) |
"""
const CTX_CPG = UInt8(0)
const CTX_CHG = UInt8(1)
const CTX_CHH = UInt8(2)
const CTX_UNKNOWN = UInt8(3)

"""Context label for each context code, in code order."""
const CONTEXT_LABELS = (:CpG, :CHG, :CHH, :unknown)

"""Narrow `value` to the two bits a payload's context or strand field holds."""
@inline two_bit_code(value::Integer) = UInt8(value) & 0x03

"""
A single genomic cytosine's aggregated methylation calls, packed into 8 bytes.

- `pos`: 1-based genomic coordinate. Positions are partitioned by scaffold (see
  [`MethylationData`](@ref)), so 32 bits is ample for any single sequence.
- `payload`: a bit-packed methylation *level* and read depth, plus metadata:

```
|-Rsvd-|--St-|--Cx-|--percent--|------------depth------------|
| 31-28|27-26|25-24|   23-16   |             15-0            |
```

- Bits  0-15 : total read depth (0-65535, saturating — see [`MAX_COUNT`](@ref))
- Bits 16-23 : methylation percentage, linearly mapped from `[0, 100]` onto
  `[0, 255]` (see [`encode_percent`](@ref))
- Bits 24-25 : cytosine context (see [`CTX_CPG`](@ref))
- Bits 26-27 : strand, as the package-wide 2-bit strand codes (see `BitCodes`)
- Bits 28-31 : reserved for future quality/SNP flags

The methylated and unmethylated counts are *reconstructed* rather than stored
(see [`get_meth`](@ref)): the reconstruction is exact up to a depth of 255 and
off by at most a read or two beyond that, while the depth stays exact and the
level stays within half a quantization step (0.196 percentage points).

Read IDs from the source alignment are deliberately discarded.
"""
struct AggregatedCall
    pos::UInt32
    payload::UInt32
end

"""
Per-scaffold aggregated methylation calls.

Each value is a `StructArray{AggregatedCall}` sorted ascending by [`sort_key`](@ref)
— position, then context, then strand — which is what lets
[`find_calls_in_range`](@ref) binary search it. Scaffolds are kept separate so a
single sequence's calls stay one contiguous, memory-mappable block.
"""
struct MethylationData
    scaffolds::Dict{String,StructArray{AggregatedCall}}
end

MethylationData() = MethylationData(Dict{String,StructArray{AggregatedCall}}())

#= Encoding / decoding =#

"""
    encode_percent(percent)

Map a methylation percentage in `[0, 100]` onto the payload's 8-bit field,
scaling linearly onto `[0, MAX_PERCENT_CODE]` and rounding half up. Values
outside `[0, 100]` clamp to the ends of the range, and `NaN` — the level
[`meth_percent`](@ref) reports for an uncovered site — encodes as `0`.
"""
@inline function encode_percent(percent::Real)
    (isnan(percent) || percent <= 0) && return UInt8(0)
    percent >= 100 && return MAX_PERCENT_CODE
    return floor(UInt8, percent * MAX_PERCENT_CODE / 100 + 0.5)
end

"""
    decode_percent(code)

Inverse of [`encode_percent`](@ref): the percentage a stored 8-bit code stands
for, as a `Float64` in `[0, 100]`. Codes are `100 / 255` ≈ 0.392 percentage
points apart, so a decoded level is within half of that of the true one.
"""
@inline decode_percent(code::Integer) = (code % Int) * 100 / MAX_PERCENT_CODE

"""
    percent_code(meth_count, unmeth_count)

Encode a pair of non-negative counts as a percentage code in integer arithmetic:
`round(MAX_PERCENT_CODE * meth / (meth + unmeth))`, half up, matching
[`encode_percent`](@ref).
"""
@inline function percent_code(meth_count::Int64, unmeth_count::Int64)
    meth_count == 0 && unmeth_count == 0 && return UInt8(0)
    # The exact form needs `510 * meth` and `2 * (meth + unmeth)` to fit in an
    # Int64; counts that large cannot come from a real library, so fall back to
    # floating point rather than risk an overflow.
    (meth_count > typemax(Int32) || unmeth_count > typemax(Int32)) && return encode_percent(
        100 * (Float64(meth_count) / (Float64(meth_count) + Float64(unmeth_count))),
    )
    total = meth_count + unmeth_count
    return UInt8(div(2 * Int64(MAX_PERCENT_CODE) * meth_count + total, 2 * total))
end

"""
    count_from_percent(code, depth)

Inverse of [`percent_code`](@ref) given a depth: `round(depth * code / 255)`,
half up. Exact for any depth ≤ 255 — see [`AggregatedCall`](@ref).
"""
@inline function count_from_percent(code::UInt8, depth::UInt32)
    scale = Int64(MAX_PERCENT_CODE)
    return UInt32(div(2 * Int64(depth) * Int64(code) + scale, 2 * scale))
end

"""
    saturating_sum(meth_count, unmeth_count)

Total depth without overflowing. The depth field tops out at
[`MAX_COUNT`](@ref), so a sum where either count alone reaches it saturates
regardless.
"""
@inline saturating_sum(meth_count::Int64, unmeth_count::Int64) =
    (meth_count >= Int64(MAX_COUNT) || unmeth_count >= Int64(MAX_COUNT)) ?
    Int64(MAX_COUNT) : meth_count + unmeth_count

"""
    pack_fields(code, depth, context, strand)

Lay an already-encoded percentage code, a depth and the metadata into the 32-bit
payload; the shared tail of [`pack_payload`](@ref) and
[`pack_percent_payload`](@ref). `depth` saturates at [`MAX_COUNT`](@ref).
"""
@inline function pack_fields(code::UInt8, depth::Integer, context::Integer, strand::Integer)
    payload = set_field(
        UInt32(0),
        clamp_field(depth, UInt32, DEPTH_WIDTH),
        DEPTH_SHIFT,
        DEPTH_WIDTH,
    )
    payload = set_field(payload, code, PERCENT_SHIFT, PERCENT_WIDTH)
    payload = set_field(payload, context, CONTEXT_SHIFT, CONTEXT_WIDTH)
    return set_field(payload, strand, STRAND_SHIFT, STRAND_FIELD_WIDTH)
end

"""
    pack_payload(meth, unmeth, context, strand)

Pack read counts and metadata into an [`AggregatedCall`](@ref) payload. The
counts are stored as a total depth plus a methylation percentage rather than
verbatim (see [`AggregatedCall`](@ref)); the percentage is computed from the
counts as given, so it stays accurate even when the depth saturates. Negative
counts are floored at zero.
"""
function pack_payload(
    meth::Integer,
    unmeth::Integer,
    context::Integer = CTX_CPG,
    strand::Integer = STRAND_NA,
)
    meth_count = max(Int64(meth), Int64(0))
    unmeth_count = max(Int64(unmeth), Int64(0))
    return pack_fields(
        percent_code(meth_count, unmeth_count),
        saturating_sum(meth_count, unmeth_count),
        context,
        strand,
    )
end

"""
    pack_percent_payload(percent, depth, context, strand)

Pack a methylation *level* and a read depth into an [`AggregatedCall`](@ref)
payload directly — the constructor to use when the level is what you have (a
model fit, a smoothed estimate, a percentage column) rather than a pair of
counts. `percent` is clamped to `[0, 100]` and `depth` to `[0, MAX_COUNT]`.
"""
pack_percent_payload(
    percent::Real,
    depth::Integer,
    context::Integer = CTX_CPG,
    strand::Integer = STRAND_NA,
) = pack_fields(encode_percent(percent), depth, context, strand)

"""
    AggregatedCall(pos, meth, unmeth, context, strand)

Build a call at `pos` from unpacked counts and metadata (see
[`pack_payload`](@ref)).
"""
AggregatedCall(
    pos::Integer,
    meth::Integer,
    unmeth::Integer,
    context::Integer = CTX_CPG,
    strand::Integer = STRAND_NA,
) = AggregatedCall(UInt32(pos), pack_payload(meth, unmeth, context, strand))

"""
Total read depth at a site, exact up to [`MAX_COUNT`](@ref). A depth of exactly
`MAX_COUNT` means the field saturated and the true depth is only known to be at
least that.
"""
@inline get_depth(payload::UInt32) = get_field(payload, DEPTH_SHIFT, DEPTH_WIDTH)

"""Raw 8-bit methylation-percentage code carried by a payload."""
@inline get_percent_code(payload::UInt32) =
    UInt8(get_field(payload, PERCENT_SHIFT, PERCENT_WIDTH))

"""
Methylation level at a site as a percentage in `[0, 100]`, or `NaN` when the
site has no coverage at all. This is the value the payload actually stores —
prefer it (or [`meth_fraction`](@ref)) over the reconstructed counts.
"""
@inline function meth_percent(payload::UInt32)
    get_depth(payload) == 0 && return NaN
    return decode_percent(get_percent_code(payload))
end

"""
Methylated read count at a site.

**Reconstructed, not stored**: the payload holds a depth and a percentage, and
this is `round(depth * percent)`. It is exact for any site whose depth is 255 or
less, and off by at most a read or two beyond that. `get_meth(p) +
get_unmeth(p) == get_depth(p)` holds for every payload, since
[`get_unmeth`](@ref) is defined as the remainder.
"""
@inline get_meth(payload::UInt32) =
    count_from_percent(get_percent_code(payload), get_depth(payload))

"""
Unmethylated read count at a site: the depth less the reconstructed methylated
count (see [`get_meth`](@ref)).
"""
@inline get_unmeth(payload::UInt32) = get_depth(payload) - get_meth(payload)

"""Cytosine context code carried by a payload (see [`CTX_CPG`](@ref))."""
@inline get_context(payload::UInt32) =
    UInt8(get_field(payload, CONTEXT_SHIFT, CONTEXT_WIDTH))

"""Context of a payload as a label: `:CpG`, `:CHG`, `:CHH` or `:unknown`."""
@inline context_label(payload::UInt32) = CONTEXT_LABELS[get_context(payload)+1]

"""Canonical 2-bit strand code carried by a payload (see `BitCodes`)."""
@inline get_strand_code(payload::UInt32) =
    UInt8(get_field(payload, STRAND_SHIFT, STRAND_FIELD_WIDTH))

"""Strand of a payload, as a `GenomicFeatures.Strand`."""
@inline BitCodes.get_strand(payload::UInt32) = decode_strand(get_strand_code(payload))

"""Whether a payload's calls came from the forward (`+`) strand."""
@inline is_forward(payload::UInt32) = get_strand_code(payload) == STRAND_FWD

"""
Fraction of reads at a site that were methylated, in `[0, 1]`, or `NaN` when the
site has no coverage at all. The stored level rescaled — see
[`meth_percent`](@ref).
"""
@inline function meth_fraction(payload::UInt32)
    get_depth(payload) == 0 && return NaN
    return get_percent_code(payload) / MAX_PERCENT_CODE
end

# Every accessor also works directly on a call.
for accessor in (
    :get_meth,
    :get_unmeth,
    :get_depth,
    :get_percent_code,
    :get_context,
    :context_label,
    :get_strand_code,
    :is_forward,
    :meth_percent,
    :meth_fraction,
)
    @eval @inline $accessor(call::AggregatedCall) = $accessor(call.payload)
end

@inline BitCodes.get_strand(call::AggregatedCall) = get_strand(call.payload)

"""Ordering key keeping call arrays in (position, context, strand) order."""
@inline sort_key(call::AggregatedCall) =
    (UInt64(call.pos) << 32) | (UInt64(get_context(call)) << 2) |
    UInt64(get_strand_code(call))

#= StructArray construction and search =#

"""
    alloc_calls(n_sites = 0)

Allocate an uninitialized `StructArray{AggregatedCall}` sized for `n_sites`. The
two 4-byte columns are allocated separately (struct of arrays), so `calls.pos`
is a contiguous `Vector{UInt32}` that `searchsorted` can scan without touching
the payloads.
"""
alloc_calls(n_sites::Integer = 0) = StructArray{AggregatedCall}((
    pos = Vector{UInt32}(undef, n_sites),
    payload = Vector{UInt32}(undef, n_sites),
))

"""
    aggregated_calls(pos, payload)

Wrap existing position and payload columns into a `StructArray{AggregatedCall}`
without copying. Both columns must have the same length; `pos` is expected to be
sorted ascending for [`find_calls_in_range`](@ref) to work.
"""
function aggregated_calls(pos::AbstractVector{UInt32}, payload::AbstractVector{UInt32})
    length(pos) == length(payload) || throw(
        DimensionMismatch(
            "pos ($(length(pos))) and payload ($(length(payload))) must have equal length.",
        ),
    )
    return StructArray{AggregatedCall}((pos = pos, payload = payload))
end

"""
    find_calls_in_range(calls, start_pos, stop_pos)

Return the calls whose position falls inside the closed interval
`[start_pos, stop_pos]`, as a **view** into `calls`: another
`StructArray{AggregatedCall}` whose columns are views of the originals, so
nothing is copied and the result can itself be range-queried.

Binary searches the contiguous `pos` column, so a lookup costs `O(log n)`.
Returns an empty view when the range holds no calls or when
`stop_pos < start_pos`.
"""
function find_calls_in_range(
    calls::StructArray{AggregatedCall},
    start_pos::Integer,
    stop_pos::Integer,
)
    stop_pos < start_pos && return view(calls, 1:0)
    positions = calls.pos
    return view(
        calls,
        searchsortedfirst(positions, start_pos):searchsortedlast(positions, stop_pos),
    )
end

"""
    find_calls_in_range(data, scaffold, start_pos, stop_pos)

Range query against one scaffold of a [`MethylationData`](@ref). Returns
`nothing` when `scaffold` isn't present.
"""
function find_calls_in_range(
    data::MethylationData,
    scaffold::AbstractString,
    start_pos::Integer,
    stop_pos::Integer,
)
    calls = get(data.scaffolds, String(scaffold), nothing)
    calls === nothing && return nothing
    return find_calls_in_range(calls, start_pos, stop_pos)
end

Base.getindex(data::MethylationData, scaffold::AbstractString) =
    data.scaffolds[String(scaffold)]
Base.haskey(data::MethylationData, scaffold::AbstractString) =
    haskey(data.scaffolds, String(scaffold))
Base.keys(data::MethylationData) = keys(data.scaffolds)
Base.length(data::MethylationData) = length(data.scaffolds)

"""Total number of aggregated sites across every scaffold."""
n_sites(data::MethylationData) = sum(length, values(data.scaffolds); init = 0)

#= Bismark parsing =#

# Byte values the line/field scanners look for.
const NEWLINE = UInt8('\n')
const CARRIAGE_RETURN = UInt8('\r')
const TAB = UInt8('\t')
const ZERO_DIGIT = UInt8('0')

"""
    parse_uint32(buffer, first_index, last_index)

Parse a non-negative *decimal* integer out of `buffer[first_index:last_index]`.
Returns `nothing` for an empty range, a non-digit byte, or a value exceeding
`typemax(UInt32)`.

Stricter than `tryparse(UInt32, ...)`, which also accepts a sign and surrounding
whitespace; neither occurs in Bismark output, so lines carrying them are dropped
as malformed.
"""
@inline function parse_uint32(
    buffer::AbstractVector{UInt8},
    first_index::Int,
    last_index::Int,
)
    first_index > last_index && return nothing
    value = UInt64(0)
    @inbounds for index = first_index:last_index
        digit = buffer[index] - ZERO_DIGIT
        digit > 0x09 && return nothing
        value = value * 10 + digit
        value > UInt64(typemax(UInt32)) && return nothing
    end
    return value % UInt32
end

"""
    decode_call(letter)

Split a Bismark methylation-call letter into `(context, is_methylated)`;
uppercase is methylated. Returns `nothing` for an unrecognised letter.
"""
@inline function decode_call(letter::Char)
    letter == 'Z' && return (CTX_CPG, true)
    letter == 'z' && return (CTX_CPG, false)
    letter == 'X' && return (CTX_CHG, true)
    letter == 'x' && return (CTX_CHG, false)
    letter == 'H' && return (CTX_CHH, true)
    letter == 'h' && return (CTX_CHH, false)
    letter == 'U' && return (CTX_UNKNOWN, true)
    letter == 'u' && return (CTX_UNKNOWN, false)
    return nothing
end

"""
    infer_strand(path)

Infer which strand a Bismark methylation-extractor file reports on from its
name: `OT`/`CTOT` files report cytosines on the forward strand, `OB`/`CTOB`
files on the reverse strand. Returns [`STRAND_NA`](@ref) when the name matches
neither.

The strand is *only* recorded in the file name — the file's own second column is
the methylation state (`+` = methylated, `-` = unmethylated), not a strand.
"""
function infer_strand(path::AbstractString)
    name = basename(String(path))
    (occursin("CTOT", name) || occursin("OT_", name) || occursin("_OT", name)) &&
        return STRAND_FWD
    (occursin("CTOB", name) || occursin("OB_", name) || occursin("_OB", name)) &&
        return STRAND_REV
    return STRAND_NA
end

"""
    parse_bismark_line(line)

Parse one Bismark methylation-extractor record into
`(scaffold, position, context, is_methylated)`, with the scaffold as a
`SubString` so nothing is allocated per record.

Expected layout (tab separated):

```
<read id>	<methylation state>	<chromosome>	<position>	<call letter>
```

Returns `nothing` for the version header, blank lines, and any line without the
five expected fields; those are skipped rather than aborting a load.
"""
@inline function parse_bismark_line(line::AbstractString)
    isempty(line) && return nothing

    tab_after_read_id = findfirst('\t', line)
    tab_after_read_id === nothing && return nothing
    tab_after_state = findnext('\t', line, tab_after_read_id + 1)
    tab_after_state === nothing && return nothing
    tab_after_scaffold = findnext('\t', line, tab_after_state + 1)
    tab_after_scaffold === nothing && return nothing
    tab_after_position = findnext('\t', line, tab_after_scaffold + 1)
    tab_after_position === nothing && return nothing

    tab_after_scaffold > tab_after_state + 1 || return nothing
    tab_after_position > tab_after_scaffold + 1 || return nothing
    lastindex(line) >= tab_after_position + 1 || return nothing

    position = parse_uint32(codeunits(line), tab_after_scaffold + 1, tab_after_position - 1)
    position === nothing && return nothing

    call = decode_call(line[tab_after_position+1])
    call === nothing && return nothing

    scaffold = SubString(line, tab_after_state + 1, tab_after_scaffold - 1)
    return (scaffold, position, call[1], call[2])
end

"""
    site_key(position, context, strand, is_methylated)

Sort/group key for one raw methylation-extractor call:

```
|--- pos (bits 63-32) ---|--- unused (31-5) ---|-cx (4-3)-|-st (2-1)-|-meth (0)-|
```

Sorting these orders records by position, and every record for one site forms a
contiguous run sharing `key >>> 1`, which [`aggregate_keys!`](@ref) collapses in
a single scan with no hashing.
"""
@inline site_key(position::UInt32, context::UInt8, strand::UInt8, is_methylated::Bool) =
    (UInt64(position) << 32) | (UInt64(context) << 3) | (UInt64(strand) << 1) |
    UInt64(is_methylated)

"""Position held by a [`site_key`](@ref)."""
@inline key_position(key::UInt64) = UInt32(key >>> 32)

"""Context code held by a [`site_key`](@ref)."""
@inline key_context(key::UInt64) = UInt8((key >>> 3) & 0x03)

"""Strand code held by a [`site_key`](@ref)."""
@inline key_strand(key::UInt64) = UInt8((key >>> 1) & 0x03)

"""Whether a [`site_key`](@ref) records a methylated call."""
@inline key_is_methylated(key::UInt64) = (key & 0x01) == 0x01

"""
    aggregate_keys!(keys)

Sort one scaffold's [`site_key`](@ref)s in place and collapse each
(position, context, strand) run into an [`AggregatedCall`](@ref).
"""
function aggregate_keys!(keys::Vector{UInt64})
    sort!(keys)

    positions = UInt32[]
    payloads = UInt32[]
    n_keys = length(keys)

    index = 1
    while index <= n_keys
        site = keys[index] >>> 1
        meth_count = 0
        unmeth_count = 0

        # Walk the run of raw calls belonging to this (pos, context, strand).
        while index <= n_keys && (keys[index] >>> 1) == site
            key_is_methylated(keys[index]) ? (meth_count += 1) : (unmeth_count += 1)
            index += 1
        end

        key = site << 1
        push!(positions, key_position(key))
        push!(
            payloads,
            pack_payload(meth_count, unmeth_count, key_context(key), key_strand(key)),
        )
    end

    return aggregated_calls(positions, payloads)
end

"""
    load_bismark(io; strand = STRAND_NA)

Read Bismark methylation-extractor records from `io` and aggregate them into
per-site counts (see [`load_bismark(::AbstractString)`](@ref) for the file-based
method, which also infers `strand` from the file name).
"""
function load_bismark(io::IO; strand::Integer = STRAND_NA)
    strand_bits = two_bit_code(strand)
    keys_by_scaffold = Dict{String,Vector{UInt64}}()

    for line in eachline(io)
        record = parse_bismark_line(line)
        record === nothing && continue
        scaffold, position, context, is_methylated = record

        # A SubString key hashes as the interned String, so a scaffold name is
        # materialised only the first time it is seen.
        bucket = get(keys_by_scaffold, scaffold, nothing)
        if bucket === nothing
            bucket = UInt64[]
            keys_by_scaffold[String(scaffold)] = bucket
        end
        push!(bucket, site_key(position, context, strand_bits, is_methylated))
    end

    return MethylationData(
        Dict{String,StructArray{AggregatedCall}}(
            scaffold => aggregate_keys!(bucket) for (scaffold, bucket) in keys_by_scaffold
        ),
    )
end

"""
    load_bismark(path; strand = infer_strand(path))

Load a Bismark methylation-extractor file (optionally gzipped) into a
[`MethylationData`](@ref) of per-site counts.

Each input line is one *read's* call at one cytosine:

```
<read id>	<methylation state>	<chromosome>	<position>	<call letter>
```

Read IDs are discarded and calls at the same site are summed, so the result is
one 8-byte [`AggregatedCall`](@ref) per (position, context, strand). The call
letter supplies the context and the methylation state (`Z`/`z` = CpG,
`X`/`x` = CHG, `H`/`h` = CHH, `U`/`u` = unknown; uppercase means methylated),
while `strand` comes from the file name — see [`infer_strand`](@ref). Pass
`strand` explicitly for files whose names don't follow the convention.

The version header, blank lines, and any malformed line are skipped.
"""
load_bismark(path::AbstractString; strand::Integer = infer_strand(path)) =
    open_maybe_gzip(io -> load_bismark(io; strand = strand), path)

"""
    open_maybe_gzip(f, path)

Open `path` and hand the stream to `f`, decompressing when the name ends in
`.gz`.
"""
function open_maybe_gzip(f, path::AbstractString)
    return open(path) do file
        io = endswith(path, ".gz") ? GzipDecompressorStream(file) : file
        try
            f(io)
        finally
            io === file || close(io)
        end
    end
end

"""
    load_bismark(paths; strand = infer_strand)

Load several Bismark files (for example the `OT` and `OB` halves of one sample)
and merge them into a single [`MethylationData`](@ref). Calls for the same
position, context *and* strand are summed; the same position on opposite strands
stays as two separate entries.

`strand` may be a fixed code applied to every file, or a function mapping a path
to a code (the default infers it from each file's name).
"""
function load_bismark(paths::AbstractVector{<:AbstractString}; strand = infer_strand)
    isempty(paths) && return MethylationData()
    return merge_calls(
        load_bismark(path; strand = strand isa Function ? strand(path) : strand) for
        path in paths
    )
end

#= Bismark coverage (.cov) parsing =#

"""
    find_tab(buffer, from_index, line_stop)

Index of the next tab in `buffer[from_index:line_stop]`, or `nothing` when the
range holds none.
"""
@inline function find_tab(buffer::Vector{UInt8}, from_index::Int, line_stop::Int)
    @inbounds for index = from_index:line_stop
        buffer[index] == TAB && return index
    end
    return nothing
end

"""
    skip_fields(buffer, from_index, line_stop, n_fields)

Index of the tab ending the `n_fields`-th field starting at `from_index`, or
`nothing` when the line runs out of fields first.
"""
@inline function skip_fields(
    buffer::Vector{UInt8},
    from_index::Int,
    line_stop::Int,
    n_fields::Int,
)
    tab = 0
    index = from_index
    for _ = 1:n_fields
        found = find_tab(buffer, index, line_stop)
        found === nothing && return nothing
        tab = found
        index = found + 1
    end
    return tab
end

"""
    parse_cov_record(buffer, line_start, line_stop)

Parse one Bismark coverage record out of `buffer[line_start:line_stop]` (a line
with its newline already stripped) into
`(scaffold_stop, position, meth_count, unmeth_count)`. The scaffold name is
`buffer[line_start:scaffold_stop]`, left as a byte range so
[`parse_cov_chunk`](@ref) can match it without allocating.

Expected layout (tab separated):

```
<chromosome>	<start>	<end>	<methylation %>	<count meth>	<count unmeth>
```

`start` is 1-based and equals `end` for a single cytosine, so it is the
position. The percentage is derivable from the counts and is skipped, as is any
column past the sixth. Returns `nothing` for blank lines and anything without
the six expected fields; those are skipped rather than aborting a load.
"""
@inline function parse_cov_record(buffer::Vector{UInt8}, line_start::Int, line_stop::Int)
    line_start > line_stop && return nothing

    tab_after_scaffold = find_tab(buffer, line_start, line_stop)
    tab_after_scaffold === nothing && return nothing
    tab_after_scaffold > line_start || return nothing        # empty scaffold name

    tab_after_start = find_tab(buffer, tab_after_scaffold + 1, line_stop)
    tab_after_start === nothing && return nothing

    # The end coordinate and the methylation percentage are both redundant.
    tab_after_percent = skip_fields(buffer, tab_after_start + 1, line_stop, 2)
    tab_after_percent === nothing && return nothing

    tab_after_meth = find_tab(buffer, tab_after_percent + 1, line_stop)
    tab_after_meth === nothing && return nothing
    tab_after_unmeth = find_tab(buffer, tab_after_meth + 1, line_stop)

    position = parse_uint32(buffer, tab_after_scaffold + 1, tab_after_start - 1)
    meth_count = parse_uint32(buffer, tab_after_percent + 1, tab_after_meth - 1)
    unmeth_count = parse_uint32(
        buffer,
        tab_after_meth + 1,
        tab_after_unmeth === nothing ? line_stop : tab_after_unmeth - 1,
    )
    (position === nothing || meth_count === nothing || unmeth_count === nothing) &&
        return nothing

    return (tab_after_scaffold - 1, position, meth_count, unmeth_count)
end

"""
One scaffold's coverage records in file order. `.cov` counts arrive already
aggregated, so unlike the extractor path there is nothing to tally: the three
columns accumulate and are then collapsed by [`collapse_cov!`](@ref).
"""
struct CovBuffer
    positions::Vector{UInt32}
    meth_counts::Vector{UInt32}
    unmeth_counts::Vector{UInt32}
end

CovBuffer() = CovBuffer(UInt32[], UInt32[], UInt32[])

"""Append one record to a [`CovBuffer`](@ref)."""
@inline function push_cov!(
    buffer::CovBuffer,
    position::UInt32,
    meth_count::UInt32,
    unmeth_count::UInt32,
)
    push!(buffer.positions, position)
    push!(buffer.meth_counts, meth_count)
    push!(buffer.unmeth_counts, unmeth_count)
    return buffer
end

"""
    collapse_cov!(buffer, context, strand)

Sort one scaffold's coverage records by position and sum repeated positions into
per-site [`AggregatedCall`](@ref)s, consuming `buffer`.

Both steps are usually no-ops — a normal `.cov` is sorted and holds one line per
cytosine — but neither is guaranteed (concatenated or re-ordered files repeat
sites) and [`find_calls_in_range`](@ref) needs the sort. Counts sum as `Int` and
saturate only when the payload is packed.
"""
function collapse_cov!(buffer::CovBuffer, context::UInt8, strand::UInt8)
    positions = buffer.positions
    meth_counts = buffer.meth_counts
    unmeth_counts = buffer.unmeth_counts

    if !issorted(positions)
        order = sortperm(positions)
        permute!(positions, order)
        permute!(meth_counts, order)
        permute!(unmeth_counts, order)
    end

    n_records = length(positions)
    site_positions = sizehint!(UInt32[], n_records)
    site_payloads = sizehint!(UInt32[], n_records)

    index = 1
    while index <= n_records
        site = positions[index]
        meth_total = Int(meth_counts[index])
        unmeth_total = Int(unmeth_counts[index])
        index += 1
        while index <= n_records && positions[index] == site
            meth_total += Int(meth_counts[index])
            unmeth_total += Int(unmeth_counts[index])
            index += 1
        end
        push!(site_positions, site)
        push!(site_payloads, pack_payload(meth_total, unmeth_total, context, strand))
    end

    return aggregated_calls(site_positions, site_payloads)
end

#= Parallel .cov ingestion =#

"""How much of the input one worker takes at a time."""
const COV_CHUNK_BYTES = 1 << 22   # 4 MiB

"""
    ChunkPool(chunk_bytes, capacity)

Free list of `chunk_bytes`-sized read buffers, so that reading a whole file does
not allocate (and immediately discard) hundreds of MiB of them.

`capacity` must be at least the number of buffers that can be alive at once, so
that [`recycle_chunk!`](@ref) never blocks.
"""
struct ChunkPool
    free::Channel{Vector{UInt8}}
    chunk_bytes::Int
end

ChunkPool(chunk_bytes::Int, capacity::Int) =
    ChunkPool(Channel{Vector{UInt8}}(capacity), chunk_bytes)

"""
    take_chunk!(pool)

A recycled buffer, or a fresh one when the pool is empty. Only the reader task
takes from a pool, so the emptiness check cannot race a competing take.
"""
@inline take_chunk!(pool::ChunkPool) =
    isready(pool.free) ? take!(pool.free) : Vector{UInt8}(undef, pool.chunk_bytes)

"""
    recycle_chunk!(pool, buffer)

Return `buffer` to `pool`. Undersized buffers (the copies
[`chunk_lines!`](@ref) makes for lines straddling a read) are dropped instead.
"""
@inline function recycle_chunk!(pool::ChunkPool, buffer::Vector{UInt8})
    length(buffer) == pool.chunk_bytes && put!(pool.free, buffer)
    return nothing
end

"""
One unit of parsing work: the whole lines in `bytes[start_index:stop_index]`.
`index` is the chunk's position in the file, which restores file order once the
workers are done. Bytes outside the range are stale content from an earlier use
of a pooled buffer and must never be read.
"""
struct CovChunk
    index::Int
    bytes::Vector{UInt8}
    start_index::Int
    stop_index::Int
end

"""
    name_matches(buffer, first_index, last_index, name)

Whether `buffer[first_index:last_index]` spells out `name`.
"""
@inline function name_matches(
    buffer::Vector{UInt8},
    first_index::Int,
    last_index::Int,
    name::String,
)
    (last_index - first_index + 1) == ncodeunits(name) || return false
    @inbounds for offset = 0:(last_index-first_index)
        buffer[first_index+offset] == codeunit(name, offset + 1) || return false
    end
    return true
end

"""
    chunk_lines!(channel, io, pool)

Cut `io` into [`CovChunk`](@ref)s of whole lines and `put!` them on `channel` in
file order, drawing read buffers from `pool`. Returns the number of chunks.

Each read is emitted as a byte range rather than a trimmed copy, so a chunk
costs no allocation. A line straddling two reads is the one exception: it is
copied into a small chunk of its own, emitted ahead of the read's body, and its
buffer is not recycled. A read containing no newline at all accumulates into
that copy.
"""
function chunk_lines!(channel::Channel{CovChunk}, io::IO, pool::ChunkPool)
    chunk_bytes = pool.chunk_bytes
    n_chunks = 0
    straddling = UInt8[]

    while true
        buffer = take_chunk!(pool)
        n_bytes = readbytes!(io, buffer, chunk_bytes)
        if n_bytes == 0
            recycle_chunk!(pool, buffer)
            break
        end

        last_newline = findprev(==(NEWLINE), buffer, n_bytes)
        if last_newline === nothing
            append!(straddling, view(buffer, 1:n_bytes))
            recycle_chunk!(pool, buffer)
            continue
        end

        body_start = 1
        if !isempty(straddling)
            # This read completes the line left over from the previous one.
            first_newline = findnext(==(NEWLINE), buffer, 1)::Int
            append!(straddling, view(buffer, 1:first_newline))
            n_chunks += 1
            put!(channel, CovChunk(n_chunks, straddling, 1, length(straddling)))
            straddling = UInt8[]
            body_start = first_newline + 1
        end

        # Copy the tail out before the buffer is handed off.
        append!(straddling, view(buffer, (last_newline+1):n_bytes))

        if body_start <= last_newline
            n_chunks += 1
            put!(channel, CovChunk(n_chunks, buffer, body_start, last_newline))
        else
            recycle_chunk!(pool, buffer)
        end
    end

    # Last line, if it had no trailing newline.
    if !isempty(straddling)
        n_chunks += 1
        put!(channel, CovChunk(n_chunks, straddling, 1, length(straddling)))
    end

    return n_chunks
end

"""
Shortest `.cov` line that can carry a record (`"c\\t1\\t1\\t0\\t0\\t0\\n"`), so a
chunk's byte count divided by this bounds its record count.
"""
const MIN_COV_LINE_BYTES = 16

"""
    reserve_cov!(buffer, n_records, rank)

Presize a new [`CovBuffer`](@ref) so its vectors do not re-pay the doubling ramp
in every chunk. `n_records` is the chunk's record ceiling and `rank` is how many
scaffolds the chunk has already opened: the first buffer gets the ceiling, each
later one half of the previous, so the reservations sum to under twice the
ceiling however many scaffolds a chunk holds.
"""
@inline function reserve_cov!(buffer::CovBuffer, n_records::Int, rank::Int)
    reserved = max(256, rank > 20 ? 0 : n_records >> (rank - 1))
    sizehint!(buffer.positions, reserved)
    sizehint!(buffer.meth_counts, reserved)
    sizehint!(buffer.unmeth_counts, reserved)
    return buffer
end

"""
How many of a chunk's scaffold names a lookup scans before falling back to the
`Dict` (see [`ScaffoldCache`](@ref)).
"""
const COV_NAME_SCAN_LIMIT = 64

"""
The scaffold buffers one [`parse_cov_chunk`](@ref) call has opened, plus the two
byte-matching layers that keep a lookup from allocating a `String` per record:
`last_name`/`last_buffer`, which a `.cov` grouped by scaffold hits for long
runs, and a scan of `names`, which covers a file that interleaves scaffolds.
`slots` is the `Dict` fallback for a chunk holding more scaffolds than
[`COV_NAME_SCAN_LIMIT`](@ref).
"""
mutable struct ScaffoldCache
    slots::Dict{String,Int}
    names::Vector{String}
    buffers::Vector{CovBuffer}
    last_name::String
    last_buffer::CovBuffer
    n_records::Int
end

ScaffoldCache(n_records::Int) =
    ScaffoldCache(Dict{String,Int}(), String[], CovBuffer[], "", CovBuffer(), n_records)

"""
    scaffold_buffer!(cache, buffer, name_start, name_stop)

Buffer for the scaffold named by `buffer[name_start:name_stop]`, opening (and
presizing) a new one the first time a name is seen.
"""
function scaffold_buffer!(
    cache::ScaffoldCache,
    buffer::Vector{UInt8},
    name_start::Int,
    name_stop::Int,
)
    name_matches(buffer, name_start, name_stop, cache.last_name) && return cache.last_buffer

    slot = 0
    if length(cache.names) <= COV_NAME_SCAN_LIMIT
        for candidate in eachindex(cache.names)
            if name_matches(buffer, name_start, name_stop, cache.names[candidate])
                slot = candidate
                break
            end
        end
    end

    if slot == 0
        name = String(buffer[name_start:name_stop])
        slot = get(cache.slots, name, 0)
        if slot == 0
            push!(cache.names, name)
            push!(
                cache.buffers,
                reserve_cov!(CovBuffer(), cache.n_records, length(cache.names)),
            )
            slot = length(cache.names)
            cache.slots[name] = slot
        end
    end

    cache.last_name = cache.names[slot]
    cache.last_buffer = cache.buffers[slot]
    return cache.last_buffer
end

"""
    parse_cov_chunk(chunk)

Parse the whole coverage lines in `chunk` into per-scaffold [`CovBuffer`](@ref)s,
keyed by scaffold name and holding records in the order they appear. Scaffold
names are resolved through a [`ScaffoldCache`](@ref).
"""
function parse_cov_chunk(chunk::CovChunk)
    bytes = chunk.bytes
    stop_index = chunk.stop_index
    cache = ScaffoldCache((stop_index - chunk.start_index + 1) ÷ MIN_COV_LINE_BYTES)

    line_start = chunk.start_index
    while line_start <= stop_index
        # A pooled buffer holds stale bytes past `stop_index`, so a newline
        # found beyond it is not this chunk's.
        newline = findnext(==(NEWLINE), bytes, line_start)
        at_end = newline === nothing || newline > stop_index
        line_stop = at_end ? stop_index : newline - 1
        # CRLF input.
        (line_stop >= line_start && bytes[line_stop] == CARRIAGE_RETURN) && (line_stop -= 1)

        record = parse_cov_record(bytes, line_start, line_stop)
        if record !== nothing
            scaffold_stop, position, meth_count, unmeth_count = record
            buffer = scaffold_buffer!(cache, bytes, line_start, scaffold_stop)
            push_cov!(buffer, position, meth_count, unmeth_count)
        end

        at_end && break
        line_start = newline + 1
    end

    return Dict{String,CovBuffer}(zip(cache.names, cache.buffers))
end

"""
    parse_cov_parallel(io, n_workers)

Read `io` and parse it into per-scaffold [`CovBuffer`](@ref)s.

One task reads (decompression, where there is any, is serial) while `n_workers`
tasks parse chunks off a bounded channel, each into its own scaffold buffers so
the hot loop shares nothing. Results concatenate in chunk order, which makes the
load independent of scheduling and leaves an already-sorted `.cov` sorted for
[`collapse_cov!`](@ref).

Read buffers come from a [`ChunkPool`](@ref) sized to the most that can be alive
at once: one in the reader's hand, `n_workers` queued, `n_workers` being parsed.
"""
function parse_cov_parallel(io::IO, n_workers::Int)
    n_workers = max(1, n_workers)
    channel = Channel{CovChunk}(n_workers)
    pool = ChunkPool(COV_CHUNK_BYTES, 2 * n_workers + 2)

    results = Tuple{Int,Dict{String,CovBuffer}}[]
    results_lock = ReentrantLock()

    workers = map(1:n_workers) do _
        Threads.@spawn begin
            try
                for chunk in channel
                    parsed = parse_cov_chunk(chunk)
                    recycle_chunk!(pool, chunk.bytes)
                    Base.@lock results_lock push!(results, (chunk.index, parsed))
                end
            catch err
                # Unblocks a reader waiting on a full channel.
                close(channel, err isa Exception ? err : ErrorException(string(err)))
                rethrow()
            end
        end
    end

    try
        chunk_lines!(channel, io, pool)
    finally
        close(channel)
    end
    foreach(wait, workers)

    sort!(results; by = first)

    # Total each scaffold up front, so the buffer the rest are folded into is
    # grown once.
    totals = Dict{String,Int}()
    for (_, parsed) in results, (scaffold, buffer) in parsed
        totals[scaffold] = get(totals, scaffold, 0) + length(buffer.positions)
    end

    merged = Dict{String,CovBuffer}()
    for (_, parsed) in results, (scaffold, buffer) in parsed
        existing = get(merged, scaffold, nothing)
        if existing === nothing
            merged[scaffold] = reserve_cov!(buffer, totals[scaffold], 1)
        else
            append!(existing.positions, buffer.positions)
            append!(existing.meth_counts, buffer.meth_counts)
            append!(existing.unmeth_counts, buffer.unmeth_counts)
        end
    end
    return merged
end

"""
    load_cov_stream(io, context, strand, n_workers)

Parse `io` on `n_workers` tasks and collapse the result into a
[`MethylationData`](@ref). Shared by every `load_bismark_cov` method;
`n_workers` is what lets concurrent files split the available threads.
"""
function load_cov_stream(io::IO, context::UInt8, strand::UInt8, n_workers::Int)
    buffers = parse_cov_parallel(io, n_workers)

    # Scaffolds collapse independently of one another.
    scaffold_names = collect(keys(buffers))
    collapsed = Vector{StructArray{AggregatedCall}}(undef, length(scaffold_names))
    @sync for (index, scaffold) in enumerate(scaffold_names)
        Threads.@spawn collapsed[index] = collapse_cov!(buffers[scaffold], context, strand)
    end

    return MethylationData(
        Dict{String,StructArray{AggregatedCall}}(zip(scaffold_names, collapsed)),
    )
end

"""
    load_bismark_cov(io; context = CTX_CPG, strand = STRAND_NA)

Read Bismark coverage records from `io` and pack them into per-site calls.

Parsing runs on `Threads.nthreads()` tasks over whole-line chunks of the stream,
and the per-scaffold collapse over the scaffolds, so this needs `julia -t auto`
to be worth anything. The result is the same at any thread count.
"""
load_bismark_cov(io::IO; context::Integer = CTX_CPG, strand::Integer = STRAND_NA) =
    load_cov_stream(io, two_bit_code(context), two_bit_code(strand), Threads.nthreads())

"""
    load_bismark_cov(path; context = CTX_CPG, strand = infer_strand(path))

Load a Bismark coverage file (`.cov`, optionally gzipped) into a
[`MethylationData`](@ref) of per-site counts.

This is Bismark's *aggregated* output — what `bismark2bedGraph`/
`coverage2cytosine` write, as opposed to the per-read methylation-extractor
records [`load_bismark`](@ref) reads. Each line is one cytosine:

```
<chromosome>	<start>	<end>	<methylation %>	<count methylated>	<count unmethylated>
```

Coordinates are already 1-based and `start == end` for a single cytosine, so
`start` is taken as the position and `end` ignored. The percentage column is
ignored as well, the counts being the more precise source, and those go straight
into an [`AggregatedCall`](@ref) whose depth saturates at [`MAX_COUNT`](@ref).
Repeated lines for the same position (concatenated files) are summed, and each
scaffold's calls are sorted by position.

A `.cov` file records **neither the context nor the strand**, so both are
supplied by the caller:

- `context` defaults to [`CTX_CPG`](@ref), Bismark's default coverage output
  being CpG-only. Pass the right code for a single-context `--CX` run, or
  [`CTX_UNKNOWN`](@ref) for a mixed-context one — a mixed file cannot be split
  by context after the fact.
- `strand` defaults to whatever [`infer_strand`](@ref) makes of the file name,
  which is [`STRAND_NA`](@ref) for the usual strand-pooled whole-sample file.
  Pass a code explicitly for coverage generated from one strand's extractor
  output.

Blank lines and any line without the six expected fields are skipped, as are
lines whose numeric fields are not plain runs of digits (a sign or padding
spaces count as malformed). Parsing is threaded — see
[`load_bismark_cov(::IO)`](@ref).
"""
load_bismark_cov(
    path::AbstractString;
    context::Integer = CTX_CPG,
    strand::Integer = infer_strand(path),
) = open_maybe_gzip(io -> load_bismark_cov(io; context = context, strand = strand), path)

"""
Default ceiling on how many `.cov` files [`load_bismark_cov`](@ref) reads at
once. Each concurrent file holds a full [`MethylationData`](@ref) plus its
[`ChunkPool`](@ref), so the bound is about memory, not scheduling.
"""
const MAX_CONCURRENT_COV_FILES = 4

"""
    load_bismark_cov(paths; context = CTX_CPG, strand = infer_strand,
                     max_concurrent_files = MAX_CONCURRENT_COV_FILES)

Load several Bismark coverage files and merge them into a single
[`MethylationData`](@ref), summing the counts of entries that share a position,
context and strand.

`context` and `strand` may each be a fixed code applied to every file, or a
function mapping a path to a code (`strand` defaults to inferring it from each
file's name).

Up to `max_concurrent_files` files are read at once, gated by a semaphore, and
the available threads are split between them; that overlaps the part of a load a
single file cannot parallelize, chiefly gzip decompression. Merging follows
`paths` order, so the result does not depend on which file finishes first.
"""
function load_bismark_cov(
    paths::AbstractVector{<:AbstractString};
    context = CTX_CPG,
    strand = infer_strand,
    max_concurrent_files::Integer = MAX_CONCURRENT_COV_FILES,
)
    isempty(paths) && return MethylationData()

    n_concurrent = clamp(Int(max_concurrent_files), 1, length(paths))
    workers_per_file = max(1, Threads.nthreads() ÷ n_concurrent)
    gate = Base.Semaphore(n_concurrent)

    tasks = map(paths) do path
        file_context = two_bit_code(context isa Function ? context(path) : context)
        file_strand = two_bit_code(strand isa Function ? strand(path) : strand)
        Threads.@spawn Base.acquire(gate) do
            open_maybe_gzip(path) do io
                load_cov_stream(io, file_context, file_strand, workers_per_file)
            end
        end
    end

    return merge_calls(fetch(task) for task in tasks)
end

#= Merging =#

"""
    merge_calls(datasets)

Merge several [`MethylationData`](@ref) into one, combining entries that share a
position, context and strand: their depths add (saturating at
[`MAX_COUNT`](@ref)) and their levels are pooled in proportion to those depths.

Merging goes through the reconstructed counts (see [`get_meth`](@ref)), so it is
exact for the depths those are exact at, and is not perfectly associative once a
site is deep enough for the level's quantization to bite.
"""
function merge_calls(datasets)
    scaffolds = Dict{String,StructArray{AggregatedCall}}()
    for data in datasets, (scaffold, calls) in data.scaffolds
        existing = get(scaffolds, scaffold, nothing)
        scaffolds[scaffold] = existing === nothing ? calls : merge_scaffold(existing, calls)
    end
    return MethylationData(scaffolds)
end

# Copy `calls[from_index:end]` onto the output columns of `merge_scaffold`.
@inline function _append_tail!(
    positions::Vector{UInt32},
    payloads::Vector{UInt32},
    calls::StructArray{AggregatedCall},
    from_index::Int,
)
    append!(positions, view(calls.pos, from_index:length(calls)))
    append!(payloads, view(calls.payload, from_index:length(calls)))
    return nothing
end

"""
    merge_scaffold(left, right)

Merge two [`sort_key`](@ref)-ordered call arrays for the same scaffold, summing
sites that share a position, context and strand.
"""
function merge_scaffold(
    left::StructArray{AggregatedCall},
    right::StructArray{AggregatedCall},
)
    n_left, n_right = length(left), length(right)
    positions = sizehint!(UInt32[], n_left + n_right)
    payloads = sizehint!(UInt32[], n_left + n_right)

    left_index, right_index = 1, 1
    while left_index <= n_left && right_index <= n_right
        left_call, right_call = left[left_index], right[right_index]
        left_key, right_key = sort_key(left_call), sort_key(right_call)

        if left_key < right_key
            push!(positions, left_call.pos)
            push!(payloads, left_call.payload)
            left_index += 1
        elseif right_key < left_key
            push!(positions, right_call.pos)
            push!(payloads, right_call.payload)
            right_index += 1
        else
            push!(positions, left_call.pos)
            push!(
                payloads,
                pack_payload(
                    Int(get_meth(left_call)) + Int(get_meth(right_call)),
                    Int(get_unmeth(left_call)) + Int(get_unmeth(right_call)),
                    get_context(left_call),
                    get_strand_code(left_call),
                ),
            )
            left_index += 1
            right_index += 1
        end
    end

    _append_tail!(positions, payloads, left, left_index)
    _append_tail!(positions, payloads, right, right_index)

    return aggregated_calls(positions, payloads)
end

#= Arrow I/O =#

"""
    write_methylation_arrow(filepath, calls; compress = :zstd)

Write one scaffold's calls to an Apache Arrow file as two `UInt32` columns
(`pos`, `payload`) — the same struct-of-arrays layout they have in memory.

`compress` accepts `:zstd`, `:lz4`, or `nothing` for no compression. The
tradeoff: compressed files are much smaller, but only an uncompressed one can be
memory-mapped by [`read_methylation_arrow`](@ref).
"""
function write_methylation_arrow(
    filepath::AbstractString,
    calls::StructArray{AggregatedCall};
    compress::Union{Nothing,Symbol} = :zstd,
)
    Arrow.write(filepath, calls; compress = compress)
    return filepath
end

"""
    read_methylation_arrow(filepath)

Read an Arrow file written by [`write_methylation_arrow`](@ref) back into a
`StructArray{AggregatedCall}` wrapping the Arrow columns directly rather than
copying them. An uncompressed file is memory-mapped and costs no resident memory
until its pages are touched; a compressed one is decompressed into memory.
"""
function read_methylation_arrow(filepath::AbstractString)
    table = Arrow.Table(filepath)
    columns = propertynames(table)
    (:pos in columns && :payload in columns) || throw(
        ArgumentError(
            "\"$filepath\" is not a methylation Arrow file (expected `pos` and `payload` columns, got $(collect(columns))).",
        ),
    )
    return StructArray{AggregatedCall}((pos = table.pos, payload = table.payload))
end

"""Name of the scaffold manifest written alongside a dataset's Arrow files."""
const MANIFEST_NAME = "scaffolds.tsv"

"""Reduce a scaffold name to something safe to use as a file name."""
sanitize_name(name::AbstractString) = replace(String(name), r"[^A-Za-z0-9._-]" => "_")

"""
    write_methylation(dir, data; compress = :zstd)

Write a whole [`MethylationData`](@ref) to `dir`, one Arrow file per scaffold
plus a `scaffolds.tsv` manifest mapping scaffold names to file names.

One file per scaffold keeps a range query to a single mapped file. The manifest
exists because scaffold names are not always valid file names: the file name is
a [`sanitize_name`](@ref)d version, deduplicated with a numeric suffix, and the
true name is recorded in the manifest.
"""
function write_methylation(
    dir::AbstractString,
    data::MethylationData;
    compress::Union{Nothing,Symbol} = :zstd,
)
    mkpath(dir)
    manifest = Tuple{String,String}[]
    used = Set{String}()

    for scaffold in sort!(collect(keys(data.scaffolds)))
        stem = sanitize_name(scaffold)
        # Two different scaffold names can sanitize to the same stem.
        candidate = stem
        suffix = 1
        while candidate in used
            candidate = "$(stem)_$(suffix)"
            suffix += 1
        end
        push!(used, candidate)

        filename = "$(candidate).arrow"
        write_methylation_arrow(
            joinpath(dir, filename),
            data.scaffolds[scaffold];
            compress = compress,
        )
        push!(manifest, (scaffold, filename))
    end

    open(joinpath(dir, MANIFEST_NAME), "w") do io
        for (scaffold, filename) in manifest
            println(io, scaffold, '\t', filename)
        end
    end

    return dir
end

"""
    read_methylation(dir)

Read a [`MethylationData`](@ref) written by [`write_methylation`](@ref) back
from `dir`, memory-mapping each scaffold's Arrow file (see
[`read_methylation_arrow`](@ref) for the caveat about compressed files).
"""
function read_methylation(dir::AbstractString)
    manifest_path = joinpath(dir, MANIFEST_NAME)
    isfile(manifest_path) ||
        throw(ArgumentError("No methylation manifest found at \"$manifest_path\"."))

    scaffolds = Dict{String,StructArray{AggregatedCall}}()
    for line in eachline(manifest_path)
        tab = findfirst('\t', line)
        tab === nothing && continue
        scaffolds[String(SubString(line, 1, tab - 1))] =
            read_methylation_arrow(joinpath(dir, String(SubString(line, tab + 1))))
    end

    return MethylationData(scaffolds)
end

#= Base.show overloads =#

function Base.show(io::IO, call::AggregatedCall)
    depth = Int(get_depth(call))
    # Show what is actually stored — a level over a depth — rather than the
    # reconstructed counts.
    level =
        depth == 0 ? "no coverage" :
        "$(round(meth_percent(call), digits = 1))% meth over $(depth) read$(depth == 1 ? "" : "s")"
    print(
        io,
        "AggregatedCall(",
        Int(call.pos),
        ", ",
        level,
        ", ",
        context_label(call),
        ", ",
        strand_char(get_strand_code(call)),
        ")",
    )
end

function Base.show(io::IO, data::MethylationData)
    n_scaffolds = length(data.scaffolds)
    sites = n_sites(data)
    print(
        io,
        "MethylationData($(n_scaffolds) scaffold$(n_scaffolds == 1 ? "" : "s"), $(sites) site$(sites == 1 ? "" : "s"))",
    )
end

export AggregatedCall,
    MethylationData,
    CTX_CPG,
    CTX_CHG,
    CTX_CHH,
    CTX_UNKNOWN,
    MAX_COUNT,
    MAX_PERCENT_CODE,
    aggregated_calls,
    alloc_calls,
    context_label,
    decode_percent,
    encode_percent,
    find_calls_in_range,
    get_context,
    get_depth,
    get_meth,
    get_unmeth,
    infer_strand,
    is_forward,
    load_bismark,
    load_bismark_cov,
    merge_calls,
    meth_fraction,
    meth_percent,
    n_sites,
    pack_payload,
    pack_percent_payload,
    read_methylation,
    read_methylation_arrow,
    write_methylation,
    write_methylation_arrow

end
