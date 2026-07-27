module Methylation

using Arrow
using CodecZlib
using StructArrays

using ..BitCodes

#= Bit layout =#

# Field layout of the 32-bit call payload (see `AggregatedCall`). Offsets are
# 0-based bit positions.
const METH_SHIFT = 0
const METH_WIDTH = 12
const UNMETH_SHIFT = 12
const UNMETH_WIDTH = 12
const CONTEXT_SHIFT = 24
const CONTEXT_WIDTH = 2
const STRAND_SHIFT = 26
const STRAND_FIELD_WIDTH = STRAND_WIDTH   # 2 bits; bits 28-31 stay reserved

"""
Largest read count a 12-bit coverage field can hold. Counts above this saturate
here rather than wrapping into a neighbouring field.
"""
const MAX_COUNT = UInt32(4095)

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
- `payload`: bit-packed counts and metadata:

```
|-Rsvd-|-St-|-Cx-|------unmeth_count------|-------meth_count-------|
| 31-28| 27 | 25 |          23-12         |          11-0          |
|  26  | 24 |
```

- Bits  0-11 : methylated read count (0-4095, saturating)
- Bits 12-23 : unmethylated read count (0-4095, saturating)
- Bits 24-25 : cytosine context (see [`CTX_CPG`](@ref))
- Bits 26-27 : strand, as the package-wide 2-bit strand codes (see `BitCodes`)
- Bits 28-31 : reserved for future quality/SNP flags

Read IDs from the source alignment are deliberately discarded: only positional
counts are kept, which is what makes the fixed 8-byte-per-site footprint (and
the columnar Arrow layout) possible.
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
    pack_payload(meth, unmeth, context, strand)

Pack read counts and metadata into an [`AggregatedCall`](@ref) payload. Counts
above [`MAX_COUNT`](@ref) are clamped (they saturate rather than wrapping into
the neighbouring field), as are negative counts (to zero).
"""
function pack_payload(
    meth::Integer,
    unmeth::Integer,
    context::Integer = CTX_CPG,
    strand::Integer = STRAND_NA,
)
    payload = UInt32(0)
    payload =
        set_field(payload, clamp_field(meth, UInt32, METH_WIDTH), METH_SHIFT, METH_WIDTH)
    payload = set_field(
        payload,
        clamp_field(unmeth, UInt32, UNMETH_WIDTH),
        UNMETH_SHIFT,
        UNMETH_WIDTH,
    )
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

"""Methylated read count carried by a payload."""
@inline get_meth(payload::UInt32) = get_field(payload, METH_SHIFT, METH_WIDTH)

"""Unmethylated read count carried by a payload."""
@inline get_unmeth(payload::UInt32) = get_field(payload, UNMETH_SHIFT, UNMETH_WIDTH)

"""
Total read depth at a site: methylated plus unmethylated. Each count saturates
at [`MAX_COUNT`](@ref) independently, so a depth of `2 * MAX_COUNT` means both
fields are saturated and the true depth is only known to be at least that.
"""
@inline get_depth(payload::UInt32) = get_meth(payload) + get_unmeth(payload)

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
Fraction of reads at a site that were methylated, or `NaN` when the site has no
coverage at all.
"""
@inline function meth_fraction(payload::UInt32)
    depth = get_depth(payload)
    return depth == 0 ? NaN : get_meth(payload) / depth
end

# Every accessor also works directly on a call.
for f in (
    :get_meth,
    :get_unmeth,
    :get_depth,
    :get_context,
    :context_label,
    :get_strand_code,
    :is_forward,
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
    return open(path) do fh
        io = endswith(path, ".gz") ? GzipDecompressorStream(fh) : fh
        try
            load_bismark(io; strand = strand)
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

"""
    merge_calls(datasets)

Merge several [`MethylationData`](@ref) into one, summing the counts of entries
that share a position, context and strand. Counts still saturate at
[`MAX_COUNT`](@ref).
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
    print(
        io,
        "AggregatedCall(",
        Int(call.pos),
        ", ",
        Int(get_meth(call)),
        "/",
        Int(get_depth(call)),
        " meth, ",
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
    aggregated_calls,
    alloc_calls,
    context_label,
    find_calls_in_range,
    get_context,
    get_depth,
    get_meth,
    get_unmeth,
    infer_strand,
    is_forward,
    load_bismark,
    merge_calls,
    meth_fraction,
    n_sites,
    pack_payload,
    read_methylation,
    read_methylation_arrow,
    write_methylation,
    write_methylation_arrow

end
