"""
Single-base methylation calls.

A submodule of [`Data`](@ref) because methylation is just another per-sample
signal sitting on top of a reference sequence. It is, however, the one such
signal that is purely *positional*: nothing here touches a `Genome`, a
`Scaffold` or any of `Data`'s interval types, and the only thing it shares with
the rest of the package is the `BitCodes` packing vocabulary. Keep it that way —
the independence is what lets a whole-genome call set stay an 8-byte-per-site,
memory-mappable block.
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
saturate here rather than wrapping into a neighbouring field.

Saturating the depth does **not** corrupt the methylation level: the percentage
is encoded from the true counts before the depth is clamped, so a site covered
by a million reads still reports its level accurately and only its depth is
recorded as "at least [`MAX_COUNT`](@ref)".
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

const CONTEXT_LABELS = (:CpG, :CHG, :CHH, :unknown)

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

The methylated and unmethylated counts are *reconstructed* rather than
stored (see [`get_meth`](@ref)), and the reconstruction is exact for every site
with a depth of 255 or less — that is, for essentially every site in a real
bisulfite library. Deeper than that the counts can be off by a read or two,
while the depth itself stays exact and the level stays within half a
quantization step (0.196 percentage points).

Read IDs from the source alignment are deliberately discarded.
"""
struct AggregatedCall
    pos::UInt32
    payload::UInt32
end

"""
Per-scaffold aggregated methylation calls.

Each value is a `StructArray{AggregatedCall}` sorted ascending by `pos` (ties —
the same position on different strands or in different contexts — are ordered by
context, then strand), which is what lets [`find_calls_in_range`](@ref) binary
search it. Scaffolds are kept separate so a single sequence's calls stay one
contiguous, memory-mappable block.
"""
struct MethylationData
    scaffolds::Dict{String,StructArray{AggregatedCall}}
end

MethylationData() = MethylationData(Dict{String,StructArray{AggregatedCall}}())

#= Encoding / decoding =#

"""
    encode_percent(percent)

Map a methylation percentage in `[0, 100]` onto the payload's 8-bit field, a
plain linear scaling onto `[0, MAX_PERCENT_CODE]` rounded (half up) to the
nearest code. Values outside `[0, 100]` clamp to the ends of the range, and
`NaN` — the level [`meth_percent`](@ref) reports for an uncovered site — encodes
as `0`, which is what an uncovered site stores anyway.
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

# Encode a pair of non-negative counts as a percentage code, staying in integer
# arithmetic: `round(MAX_PERCENT_CODE * m / (m + u))`, half up, matching
# `encode_percent`.
@inline function percent_code(m::Int64, u::Int64)
    m == 0 && u == 0 && return UInt8(0)
    # The exact form needs `510 * m` and `2 * (m + u)` to fit in an Int64.
    # Counts that large cannot come from a real library, so rather than risk an
    # overflow, fall back to floating point for them.
    (m > typemax(Int32) || u > typemax(Int32)) &&
        return encode_percent(100 * (Float64(m) / (Float64(m) + Float64(u))))
    total = m + u
    return UInt8(div(2 * Int64(MAX_PERCENT_CODE) * m + total, 2 * total))
end

# Inverse of `percent_code` given a depth: `round(depth * code / 255)`, half up.
# Exact for any depth <= 255; see the `AggregatedCall` docstring.
@inline function count_from_percent(code::UInt8, depth::UInt32)
    scale = Int64(MAX_PERCENT_CODE)
    return UInt32(div(2 * Int64(depth) * Int64(code) + scale, 2 * scale))
end

# `m + u` without ever overflowing: the depth field tops out at MAX_COUNT, so if
# either count alone already reaches it the sum saturates regardless.
@inline saturating_sum(m::Int64, u::Int64) =
    (m >= Int64(MAX_COUNT) || u >= Int64(MAX_COUNT)) ? Int64(MAX_COUNT) : m + u

"""
    pack_payload(meth, unmeth, context, strand)

Pack read counts and metadata into an [`AggregatedCall`](@ref) payload.

The counts are stored as a total depth plus a methylation percentage rather than
verbatim (see [`AggregatedCall`](@ref)), so they are recovered by
[`get_meth`](@ref)/[`get_unmeth`](@ref) exactly whenever their sum is 255 or
less. The percentage is computed from the counts as given, so it stays accurate
even when the depth itself saturates at [`MAX_COUNT`](@ref); negative counts are
floored at zero.
"""
function pack_payload(
    meth::Integer,
    unmeth::Integer,
    context::Integer = CTX_CPG,
    strand::Integer = STRAND_NA,
)
    m = max(Int64(meth), Int64(0))
    u = max(Int64(unmeth), Int64(0))
    return pack_fields(percent_code(m, u), saturating_sum(m, u), context, strand)
end

"""
    pack_percent_payload(percent, depth, context, strand)

Pack a methylation *level* and a read depth into an [`AggregatedCall`](@ref)
payload directly — the natural constructor when the level is what you have (a
model fit, a smoothed estimate, a percentage column) rather than a pair of
counts. `percent` is clamped to `[0, 100]` and `depth` to
`[0, MAX_COUNT]`.
"""
function pack_percent_payload(
    percent::Real,
    depth::Integer,
    context::Integer = CTX_CPG,
    strand::Integer = STRAND_NA,
)
    return pack_fields(encode_percent(percent), depth, context, strand)
end

# Shared tail of both packers: lay an already-encoded percentage code, a depth
# and the metadata into the 32-bit payload.
@inline function pack_fields(code::UInt8, depth::Integer, context::Integer, strand::Integer)
    payload = UInt32(0)
    payload = set_field(
        payload,
        clamp_field(depth, UInt32, DEPTH_WIDTH),
        DEPTH_SHIFT,
        DEPTH_WIDTH,
    )
    payload = set_field(payload, code, PERCENT_SHIFT, PERCENT_WIDTH)
    payload = set_field(payload, context, CONTEXT_SHIFT, CONTEXT_WIDTH)
    payload = set_field(payload, strand, STRAND_SHIFT, STRAND_FIELD_WIDTH)
    return payload
end

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
least that; the methylation level is unaffected either way.
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
for f in (
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
    @eval @inline $f(call::AggregatedCall) = $f(call.payload)
end

@inline BitCodes.get_strand(call::AggregatedCall) = get_strand(call.payload)

#= StructArray construction and search =#

"""
    alloc_calls(n = 0)

Allocate an uninitialized `StructArray{AggregatedCall}` sized for `n` sites. The
two 4-byte columns are allocated separately (struct of arrays), so `calls.pos`
is a contiguous `Vector{UInt32}` that `searchsorted` can scan without touching
the payloads.
"""
alloc_calls(n::Integer = 0) = StructArray{AggregatedCall}((
    pos = Vector{UInt32}(undef, n),
    payload = Vector{UInt32}(undef, n),
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
`StructArray{AggregatedCall}`, but one whose columns are views of the originals,
so nothing is copied (memory-mapped columns included) and the result can itself
be range-queried.

Binary searches the contiguous `pos` column, so a lookup costs `O(log n)`
regardless of how many sites the scaffold holds. Returns an empty view when the
range holds no calls or when `stop_pos < start_pos`.
"""
function find_calls_in_range(
    calls::StructArray{AggregatedCall},
    start_pos::Integer,
    stop_pos::Integer,
)
    positions = calls.pos
    stop_pos < start_pos && return view(calls, 1:0)
    lo = searchsortedfirst(positions, start_pos)
    hi = searchsortedlast(positions, stop_pos)
    return view(calls, lo:hi)
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

# Bismark methylation-call letters: uppercase = methylated, lowercase = not.
# Returns (context, is_methylated), or `nothing` for anything unrecognised.
@inline function decode_call(c::Char)
    c == 'Z' && return (CTX_CPG, true)
    c == 'z' && return (CTX_CPG, false)
    c == 'X' && return (CTX_CHG, true)
    c == 'x' && return (CTX_CHG, false)
    c == 'H' && return (CTX_CHH, true)
    c == 'h' && return (CTX_CHH, false)
    c == 'U' && return (CTX_UNKNOWN, true)
    c == 'u' && return (CTX_UNKNOWN, false)
    return nothing
end

"""
    infer_strand(path)

Infer which strand a Bismark methylation-extractor file reports on from its
name. Bismark splits its output by alignment strand, and the strand is *only*
recorded in the file name — the file's own second column is the methylation
state (`+` = methylated, `-` = unmethylated), not a strand.

`OT`/`CTOT` files report cytosines on the forward strand, `OB`/`CTOB` files
report cytosines on the reverse strand. Returns [`STRAND_NA`](@ref) when the
name matches neither.
"""
function infer_strand(path::AbstractString)
    name = basename(String(path))
    (occursin("CTOT", name) || occursin("OT_", name) || occursin("_OT", name)) &&
        return STRAND_FWD
    (occursin("CTOB", name) || occursin("OB_", name) || occursin("_OB", name)) &&
        return STRAND_REV
    return STRAND_NA
end

# Parse one Bismark record line into (chrom, pos, context, is_methylated).
#
# Expected layout (tab separated):
#   <read id> <methylation state> <chromosome> <position> <call letter>
#
# Returns `nothing` for the version header, blank lines, and any line that
# doesn't have the five expected fields, so those are skipped rather than
# aborting a load. Field boundaries are located by index so the chromosome and
# position are read as `SubString`s, without allocating per record.
@inline function parse_bismark_line(line::AbstractString)
    isempty(line) && return nothing

    t1 = findfirst('\t', line)
    t1 === nothing && return nothing
    t2 = findnext('\t', line, t1 + 1)
    t2 === nothing && return nothing
    t3 = findnext('\t', line, t2 + 1)
    t3 === nothing && return nothing
    t4 = findnext('\t', line, t3 + 1)
    t4 === nothing && return nothing

    t3 > t2 + 1 || return nothing
    t4 > t3 + 1 || return nothing
    lastindex(line) >= t4 + 1 || return nothing

    chrom = SubString(line, t2 + 1, t3 - 1)
    pos = tryparse(UInt32, SubString(line, t3 + 1, t4 - 1))
    pos === nothing && return nothing

    call = decode_call(line[t4+1])
    call === nothing && return nothing

    return (chrom, pos, call[1], call[2])
end

# Sort/group key for one raw call. Laying the fields out as
#
#   |--- pos (bits 63-32) ---|--- unused (31-5) ---|-cx (4-3)-|-st (2-1)-|-meth (0)-|
#
# means an ordinary `sort!` of the keys orders records by position first, and
# that every record for one site lands in a contiguous run whose members share
# `key >>> 1`. Aggregation is then a single linear scan, with no hashing and no
# per-site allocation.
@inline site_key(pos::UInt32, context::UInt8, strand::UInt8, meth::Bool) =
    (UInt64(pos) << 32) | (UInt64(context) << 3) | (UInt64(strand) << 1) | UInt64(meth)

@inline key_pos(key::UInt64) = UInt32(key >>> 32)
@inline key_context(key::UInt64) = UInt8((key >>> 3) & 0x03)
@inline key_strand(key::UInt64) = UInt8((key >>> 1) & 0x03)
@inline key_meth(key::UInt64) = (key & 0x01) == 0x01

# Collapse a scaffold's sorted raw keys into aggregated per-site calls.
function aggregate_keys!(keys::Vector{UInt64})
    sort!(keys)

    positions = UInt32[]
    payloads = UInt32[]

    i = 1
    n = length(keys)
    while i <= n
        site = keys[i] >>> 1
        meth = 0
        unmeth = 0

        # Walk the run of raw calls belonging to this (pos, context, strand).
        while i <= n && (keys[i] >>> 1) == site
            key_meth(keys[i]) ? (meth += 1) : (unmeth += 1)
            i += 1
        end

        key = site << 1
        push!(positions, key_pos(key))
        push!(payloads, pack_payload(meth, unmeth, key_context(key), key_strand(key)))
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
    strand_bits = UInt8(strand) & 0x03
    raw = Dict{String,Vector{UInt64}}()

    for line in eachline(io)
        record = parse_bismark_line(line)
        record === nothing && continue
        chrom, pos, context, meth = record

        # `get` with a SubString key hits the same hash as the interned String,
        # so a scaffold name is only materialised the first time it is seen.
        bucket = get(raw, chrom, nothing)
        if bucket === nothing
            bucket = UInt64[]
            raw[String(chrom)] = bucket
        end
        push!(bucket, site_key(pos, context, strand_bits, meth))
    end

    scaffolds = Dict{String,StructArray{AggregatedCall}}()
    for (chrom, bucket) in raw
        scaffolds[chrom] = aggregate_keys!(bucket)
    end
    return MethylationData(scaffolds)
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
one 8-byte [`AggregatedCall`](@ref) per (position, context, strand) rather than
one entry per read. The call letter supplies the context and the methylation
state (`Z`/`z` = CpG, `X`/`x` = CHG, `H`/`h` = CHH, `U`/`u` = unknown; uppercase
means methylated), while `strand` comes from the file name — see
[`infer_strand`](@ref) — because Bismark records it nowhere else. Pass `strand`
explicitly for files whose names don't follow the convention.

The version header, blank lines, and any malformed line are skipped.
"""
function load_bismark(path::AbstractString; strand::Integer = infer_strand(path))
    return open_maybe_gzip(path) do io
        load_bismark(io; strand = strand)
    end
end

# Open `path` for reading and hand the stream to `f`, transparently
# decompressing it when the name ends in ".gz".
function open_maybe_gzip(f, path::AbstractString)
    return open(path) do fh
        io = endswith(path, ".gz") ? GzipDecompressorStream(fh) : fh
        try
            f(io)
        finally
            io === fh || close(io)
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

# Parse one Bismark coverage line into (chrom, pos, meth, unmeth).
#
# Expected layout (tab separated):
#   <chromosome> <start> <end> <methylation %> <count meth> <count unmeth>
#
# Unlike the methylation-extractor format this is already aggregated: one line
# per cytosine, carrying counts rather than a single read's call. `start` is
# 1-based (and equals `end` for a single cytosine), so it is used as-is; the
# percentage is redundant with the two counts and is skipped. Trailing columns
# beyond the sixth are tolerated and ignored.
#
# Returns `nothing` for blank lines and anything without the six expected
# fields, so those are skipped rather than aborting a load. Field boundaries are
# located by index so the chromosome is read as a `SubString`, without
# allocating per record.
@inline function parse_cov_line(line::AbstractString)
    isempty(line) && return nothing

    t1 = findfirst('\t', line)
    t1 === nothing && return nothing
    t2 = findnext('\t', line, t1 + 1)
    t2 === nothing && return nothing
    t3 = findnext('\t', line, t2 + 1)
    t3 === nothing && return nothing
    t4 = findnext('\t', line, t3 + 1)
    t4 === nothing && return nothing
    t5 = findnext('\t', line, t4 + 1)
    t5 === nothing && return nothing

    t1 > firstindex(line) || return nothing
    lastindex(line) >= t5 + 1 || return nothing

    chrom = SubString(line, firstindex(line), t1 - 1)
    pos = tryparse(UInt32, SubString(line, t1 + 1, t2 - 1))
    pos === nothing && return nothing
    meth = tryparse(UInt32, SubString(line, t4 + 1, t5 - 1))
    meth === nothing && return nothing

    t6 = findnext('\t', line, t5 + 1)
    unmeth =
        tryparse(UInt32, SubString(line, t5 + 1, t6 === nothing ? lastindex(line) : t6 - 1))
    unmeth === nothing && return nothing

    return (chrom, pos, meth, unmeth)
end

# One scaffold's coverage records in file order. The counts arrive already
# aggregated, so unlike the extractor path there is nothing to tally — the three
# columns are just accumulated and then collapsed by `collapse_cov!`.
struct CovBuffer
    pos::Vector{UInt32}
    meth::Vector{UInt32}
    unmeth::Vector{UInt32}
end

CovBuffer() = CovBuffer(UInt32[], UInt32[], UInt32[])

@inline function push_cov!(buf::CovBuffer, pos::UInt32, meth::UInt32, unmeth::UInt32)
    push!(buf.pos, pos)
    push!(buf.meth, meth)
    push!(buf.unmeth, unmeth)
    return buf
end

# Collapse one scaffold's coverage records into position-sorted per-site calls.
#
# A `.cov` file is normally already sorted and holds one line per cytosine, so
# both the sort and the duplicate-summing scan are usually no-ops — but neither
# is guaranteed (concatenated or re-ordered files repeat sites), and
# `find_calls_in_range` needs the sort, so both are done regardless. Counts are
# summed as `Int` and only saturate once, when the payload is packed.
function collapse_cov!(buf::CovBuffer, context::UInt8, strand::UInt8)
    positions, meth, unmeth = buf.pos, buf.meth, buf.unmeth

    if !issorted(positions)
        order = sortperm(positions)
        permute!(positions, order)
        permute!(meth, order)
        permute!(unmeth, order)
    end

    n = length(positions)
    out_pos = UInt32[]
    out_payload = UInt32[]
    sizehint!(out_pos, n)
    sizehint!(out_payload, n)

    i = 1
    while i <= n
        site = positions[i]
        m = Int(meth[i])
        u = Int(unmeth[i])
        i += 1
        while i <= n && positions[i] == site
            m += Int(meth[i])
            u += Int(unmeth[i])
            i += 1
        end
        push!(out_pos, site)
        push!(out_payload, pack_payload(m, u, context, strand))
    end

    return aggregated_calls(out_pos, out_payload)
end

"""
    load_bismark_cov(io; context = CTX_CPG, strand = STRAND_NA)

Read Bismark coverage records from `io` and pack them into per-site calls.
"""
function load_bismark_cov(io::IO; context::Integer = CTX_CPG, strand::Integer = STRAND_NA)
    context_bits = UInt8(context) & 0x03
    strand_bits = UInt8(strand) & 0x03
    raw = Dict{String,CovBuffer}()

    for line in eachline(io)
        record = parse_cov_line(line)
        record === nothing && continue
        chrom, pos, meth, unmeth = record

        # `get` with a SubString key hits the same hash as the interned String,
        # so a scaffold name is only materialised the first time it is seen.
        buf = get(raw, chrom, nothing)
        if buf === nothing
            buf = CovBuffer()
            raw[String(chrom)] = buf
        end
        push_cov!(buf, pos, meth, unmeth)
    end

    scaffolds = Dict{String,StructArray{AggregatedCall}}()
    for (chrom, buf) in raw
        scaffolds[chrom] = collapse_cov!(buf, context_bits, strand_bits)
    end
    return MethylationData(scaffolds)
end

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
ignored as well — it is derivable from the two counts, which are the more
precise source — and the counts go straight into an [`AggregatedCall`](@ref),
whose depth saturates at [`MAX_COUNT`](@ref). Repeated lines for the same
position (concatenated files) are summed, and each scaffold's calls are sorted
by position.

A `.cov` file records **neither the context nor the strand**, so both are
supplied by the caller:

- `context` defaults to [`CTX_CPG`](@ref) because Bismark's default coverage
  output is CpG-only. Pass the right code for a single-context `--CX` run, or
  [`CTX_UNKNOWN`](@ref) for a mixed-context one — a mixed file cannot be split
  by context after the fact, since the information simply is not in it.
- `strand` defaults to whatever [`infer_strand`](@ref) makes of the file name,
  which is [`STRAND_NA`](@ref) for the usual whole-sample coverage file (its
  counts are pooled across strands anyway). Pass a code explicitly for coverage
  generated from one strand's extractor output.

Blank lines and any line without the six expected fields are skipped.
"""
function load_bismark_cov(
    path::AbstractString;
    context::Integer = CTX_CPG,
    strand::Integer = infer_strand(path),
)
    return open_maybe_gzip(path) do io
        load_bismark_cov(io; context = context, strand = strand)
    end
end

"""
    load_bismark_cov(paths; context = CTX_CPG, strand = infer_strand)

Load several Bismark coverage files and merge them into a single
[`MethylationData`](@ref), summing the counts of entries that share a position,
context and strand.

`context` and `strand` may each be a fixed code applied to every file, or a
function mapping a path to a code (`strand` defaults to inferring it from each
file's name).
"""
function load_bismark_cov(
    paths::AbstractVector{<:AbstractString};
    context = CTX_CPG,
    strand = infer_strand,
)
    isempty(paths) && return MethylationData()
    return merge_calls(
        load_bismark_cov(
            path;
            context = context isa Function ? context(path) : context,
            strand = strand isa Function ? strand(path) : strand,
        ) for path in paths
    )
end

"""
    merge_calls(datasets)

Merge several [`MethylationData`](@ref) into one, combining entries that share a
position, context and strand: their depths add (saturating at
[`MAX_COUNT`](@ref)) and their levels are pooled in proportion to those depths.

Merging goes through the reconstructed counts (see [`get_meth`](@ref)), so it is
exact for the depths those are exact at, and merging is not perfectly
associative once a site is deep enough for the level's quantization to bite.
"""
function merge_calls(datasets)
    scaffolds = Dict{String,StructArray{AggregatedCall}}()
    for data in datasets
        for (chrom, calls) in data.scaffolds
            existing = get(scaffolds, chrom, nothing)
            scaffolds[chrom] =
                existing === nothing ? calls : merge_scaffold(existing, calls)
        end
    end
    return MethylationData(scaffolds)
end

# Merge two sorted call arrays for the same scaffold, summing shared sites.
function merge_scaffold(a::StructArray{AggregatedCall}, b::StructArray{AggregatedCall})
    positions = UInt32[]
    payloads = UInt32[]
    sizehint!(positions, length(a) + length(b))
    sizehint!(payloads, length(a) + length(b))

    # Both inputs are sorted by (pos, context, strand), so this is a merge of two
    # ordered runs; entries that agree on all three are summed as they meet.
    i, j = 1, 1
    while i <= length(a) || j <= length(b)
        take_a = j > length(b) || (i <= length(a) && sort_key(a[i]) <= sort_key(b[j]))
        if take_a && i <= length(a) && j <= length(b) && sort_key(a[i]) == sort_key(b[j])
            call_a, call_b = a[i], b[j]
            push!(positions, call_a.pos)
            push!(
                payloads,
                pack_payload(
                    Int(get_meth(call_a)) + Int(get_meth(call_b)),
                    Int(get_unmeth(call_a)) + Int(get_unmeth(call_b)),
                    get_context(call_a),
                    get_strand_code(call_a),
                ),
            )
            i += 1
            j += 1
        elseif take_a
            push!(positions, a[i].pos)
            push!(payloads, a[i].payload)
            i += 1
        else
            push!(positions, b[j].pos)
            push!(payloads, b[j].payload)
            j += 1
        end
    end

    return aggregated_calls(positions, payloads)
end

# Ordering key used to keep merged arrays in (pos, context, strand) order.
@inline sort_key(call::AggregatedCall) =
    (UInt64(call.pos) << 32) | (UInt64(get_context(call)) << 2) |
    UInt64(get_strand_code(call))

#= Arrow I/O =#

"""
    write_methylation_arrow(filepath, calls; compress = :zstd)

Write one scaffold's calls to an Apache Arrow file as two `UInt32` columns
(`pos`, `payload`) — the same struct-of-arrays layout they have in memory.

`compress` accepts `:zstd`, `:lz4`, or `nothing` for no compression. Note the
tradeoff: compressed files are much smaller, but
[`read_methylation_arrow`](@ref) has to decompress them into memory, whereas an
uncompressed file is memory-mapped and costs no resident memory until touched.
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
`StructArray{AggregatedCall}`.

The table is memory-mapped, and the returned `StructArray` wraps the Arrow
columns directly rather than copying them, so an uncompressed file costs no
resident memory until its pages are touched. Compressed files are decompressed
into memory by Arrow on read.
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

"""
    write_methylation(dir, data; compress = :zstd)

Write a whole [`MethylationData`](@ref) to `dir`, one Arrow file per scaffold
plus a `scaffolds.tsv` manifest mapping scaffold names to file names.

Keeping each scaffold in its own file is what makes large datasets workable:
a range query only ever needs the one file to be mapped. The manifest exists
because scaffold names are not always valid file names, so the file name is a
sanitized version and the true name is recorded in the manifest.
"""
function write_methylation(
    dir::AbstractString,
    data::MethylationData;
    compress::Union{Nothing,Symbol} = :zstd,
)
    mkpath(dir)
    manifest = Tuple{String,String}[]
    used = Set{String}()

    for chrom in sort!(collect(keys(data.scaffolds)))
        stem = sanitize_name(chrom)
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
            data.scaffolds[chrom];
            compress = compress,
        )
        push!(manifest, (chrom, filename))
    end

    open(joinpath(dir, MANIFEST_NAME), "w") do io
        for (chrom, filename) in manifest
            println(io, chrom, '\t', filename)
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
        isempty(line) && continue
        tab = findfirst('\t', line)
        tab === nothing && continue
        chrom = String(SubString(line, 1, tab - 1))
        filename = String(SubString(line, tab + 1))
        scaffolds[chrom] = read_methylation_arrow(joinpath(dir, filename))
    end

    return MethylationData(scaffolds)
end

# Reduce a scaffold name to something safe to use as a file name.
sanitize_name(name::AbstractString) = replace(String(name), r"[^A-Za-z0-9._-]" => "_")

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
    n = length(data.scaffolds)
    sites = n_sites(data)
    print(
        io,
        "MethylationData($(n) scaffold$(n == 1 ? "" : "s"), $(sites) site$(sites == 1 ? "" : "s"))",
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
