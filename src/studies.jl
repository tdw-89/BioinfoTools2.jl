module Studies

using BED
using BioGenerics
using CodecZlib
using DataFrames
using Dates
using GFF3
using IntervalTrees
using SparseArrays
using ..Reference

mutable struct BedData
    scaffolds::Dict{String, IntervalMeta64}

end

mutable struct TabularData{T}
    variables::Vector{AbstractString}
    samples::Vector{Union{Nothing, Tuple{String, UInt32}}}
    table::Matrix{T}
end

const Data = Union{BedData, TabularData}

mutable struct BioSample
    sample_id::String
    tissue_type::String
    species::Species
end

mutable struct AssayMethod
    name::String
    description::String
end

mutable struct Measurement
    file_path::String
    format::String
    data::Data
end

mutable struct Assay
    id::String
    type::String
    description::String
    measurement::Measurement
    biosample::BioSample
    method::AssayMethod
end

mutable struct Study
    id::String
    title::String
    date::Date
    assays::Vector{Assay}
end

#= Methods =#

"""
NOTE: This functions assumes that the first column is a list of 
sample IDs (e.g., gene accessions or transcript accessions)
and that all other columns are of the same numeric type.
"""
function load_table(genome::Genome, data_frame::DataFrame, feature::Union{String, Symbol})
    sample_col_correct = typeof(Vector(data_frame[!, 1])) <: Array{S, 1} where S <: AbstractString
    data_cols_correct = [typeof(Vector(data_frame[!, i])) <: Array{N, 1} where N <: Number for i in 2:ncol(data_frame)] |> all
    if !(sample_col_correct && data_cols_correct)
        if !sample_col_correct
            @warn "Sample column incorrect format"
        end
        if !data_cols_correct
            @warn "Data columns incorrect format"
        end
        return nothing
    end

    # (sample_id, original 1-based row index) for every row of the table.
    sample_names = collect(zip(data_frame[!, 1], 1:nrow(data_frame)))

    # Output vector, parallel to the table rows. A slot stays `nothing` if its
    # sample ID never matches any feature's metadata ID.
    samples = Vector{Union{Nothing, Tuple{String, UInt32}}}(nothing, length(sample_names))

    # Samples still awaiting a match: sample_id => original row index. Using a
    # Dict makes each match an O(1) lookup + "pop", so the whole search is
    # O(features) and can stop early once every sample has been matched.
    remaining = Dict{String, UInt32}()
    sizehint!(remaining, length(sample_names))
    for (id, idx) in sample_names
        remaining[String(id)] = UInt32(idx)
    end

    # 1. Walk every scaffold, restricting to features of the requested type, and
    # 2. match each feature's metadata ID against the remaining samples.
    for (scaffold_name, scaffold) in genome.scaffolds
        isempty(remaining) && break

        feature_intervals = get_feature(scaffold, feature)
        (isnothing(feature_intervals) || isempty(feature_intervals)) && continue

        for interval in feature_intervals
            # The 32-bit metadata index links this interval back to its metadata.
            meta_idx = Reference.parse_index(interval.value)
            id = Reference.get_metadata_id(genome, meta_idx)
            isnothing(id) && continue

            row = get(remaining, id, nothing)
            if row !== nothing
                samples[row] = (scaffold_name, meta_idx)
                delete!(remaining, id)   # pop the matched sample
                isempty(remaining) && break
            end
        end
    end

    variables = Vector{AbstractString}(names(data_frame)[2:end])
    table = Matrix(data_frame[!, 2:end])

    return TabularData(variables, samples, table)
end

#= TabularData indexing =#

# Integer / cartesian indexing behaves exactly like indexing the underlying
# matrix (e.g. `tab[2, 2]`, `tab[1:10, :]`).
Base.getindex(t::TabularData, inds...) = getindex(t.table, inds...)

# Resolve the ID of a single sample by looking its metadata up in `genome`.
# Returns `nothing` for empty slots or features without metadata.
function sample_id(genome::Genome, sample::Union{Nothing, Tuple{String, UInt32}})
    isnothing(sample) && return nothing
    return Reference.get_metadata_id(genome, sample[2])
end

"""
Look up a single sample by its (putative) ID string. Every sample's ID is
resolved from `genome` and compared to `id`; the first match is returned as a
one-row `TabularData`. Returns `nothing` when no sample matches.

This is intentionally O(samples) per lookup (it resolves each sample's ID on
demand) rather than maintaining an ID index.
"""
function Base.getindex(t::TabularData, genome::Genome, id::AbstractString)
    for (row, sample) in enumerate(t.samples)
        if sample_id(genome, sample) == id
            return TabularData(copy(t.variables), t.samples[row:row], t.table[row:row, :])
        end
    end
    return nothing
end

"""
Look up every sample whose (putative) ID is contained in `ids`, returning a
`TabularData` sub-table (with matching `variables`, `samples` and `table`
rows). IDs with no matching sample are silently skipped; row order follows the
original table.
"""
function Base.getindex(t::TabularData, genome::Genome, ids::AbstractVector{<:AbstractString})
    wanted = Set{String}(ids)
    rows = Int[]
    for (row, sample) in enumerate(t.samples)
        sid = sample_id(genome, sample)
        if !isnothing(sid) && sid in wanted
            push!(rows, row)
        end
    end
    return TabularData(copy(t.variables), t.samples[rows], t.table[rows, :])
end

"""
Convert a `TabularData` back into a `DataFrame`. The first column, `ID`, holds
the ID that each sample's metadata index points to (resolved from `genome`;
unmatched samples become `missing`), followed by one column per entry in
`variables` carrying the corresponding `table` column.
"""
function DataFrames.DataFrame(t::TabularData, genome::Genome)
    df = DataFrame()
    df[!, :ID] = Union{Missing, String}[something(sample_id(genome, s), missing) for s in t.samples]
    for (j, variable) in enumerate(t.variables)
        df[!, string(variable)] = t.table[:, j]
    end
    return df
end

"""
Bit layout for the 64-bit metadata code attached to each BED interval.
Bits 33-40 encode strand using the same scheme as gene intervals in `Scaffold`.

|  64-57  |  56-41  | 40-33  |          32-1          |
|---------|---------|--------|------------------------|
| <NULL>  | <NULL>  | Strand |        <NULL>          |

- Bits  1-32  : reserved
- Bits 33-40  : strand (0 = unknown/unstranded, 1 = +, 2 = -)
- Bits 41-64  : reserved
"""
function pack_bed_code(strand::UInt8)
    return UInt64(strand) << 32
end

function parse_bed_strand(code::UInt64)
    return UInt8((code >> 32) & 0xFF)
end

"""
Convert a BED record's strand field to the same UInt8 encoding used by genes.
"""
function bed_record_strand(record::BED.Record)
    BED.hasstrand(record) || return UInt8(0)
    s = BED.strand(record)
    if s == GFF3.GenomicFeatures.STRAND_POS
        UInt8(1)
    elseif s == GFF3.GenomicFeatures.STRAND_NEG
        UInt8(2)
    else
        UInt8(0)
    end
end

function load_bed(file_path::String)
    scaffolds = Dict{String, IntervalMeta64}()

    open(file_path) do fh
        rdr = endswith(file_path, ".gz") ?
            BED.Reader(GzipDecompressorStream(fh)) :
            BED.Reader(fh)

        record = BED.Record()
        try
            while !eof(rdr)
                empty!(record)
                read!(rdr, record)
                # Process if filled: handles both normal reads and the last record
                # in files without a trailing newline (EOFError thrown after fill).
                if BED.isfilled(record)
                    # BED uses 0-based half-open [start, end) coordinates;
                    # convert to 1-based closed [start, end] to match gene intervals.
                    # from GFF3 files.
                    chrom     = BED.chrom(record)
                    start_pos = UInt32(BED.chromstart(record))
                    end_pos   = UInt32(BED.chromend(record))

                    strand = bed_record_strand(record)
                    code   = pack_bed_code(strand)

                    tree = get!(IntervalMeta64, scaffolds, chrom)
                    push!(tree, IntervalValue(start_pos, end_pos, code))
                end
            end
        finally
            close(rdr)
        end
    end

    return BedData(scaffolds)
end

function intersect(tree_a::IntervalMeta64, tree_b::IntervalMeta64)
    intersection = IntervalMeta64()
    itr = IntervalTrees.intersect(tree_a, tree_b)
    for overlap in itr
        a = overlap[1].first:overlap[1].last
        b = overlap[2].first:overlap[2].last
        a_x_b = Base.intersect(a,b)
        new_interval = IntervalValue(first(a_x_b), last(a_x_b), overlap[1].value)
        push!(intersection, new_interval)
    end
    return intersection
end

function intersect(
    scaffold::Scaffold,
    bed_data::BedData,
    feature::Union{AbstractString, Symbol})
    if !haskey(bed_data.scaffolds, scaffold.name)
        return nothing
    end
    feature_intervals = get_feature(scaffold, feature)
    if isnothing(feature_intervals)
        return nothing
    end
    return intersect(feature_intervals, bed_data.scaffolds[scaffold.name])
end

function intersect(
    genome::Genome,
    bed_data::BedData,
    feature::Union{AbstractString, Symbol})
    scaffolds = Dict{String, IntervalMeta64}()
    for (scaffold_name, scaffold) in genome.scaffolds
        if haskey(bed_data.scaffolds, scaffold_name)
            intersect_result = intersect(scaffold, bed_data, feature)
            if !isnothing(intersect_result)
                scaffolds[scaffold_name] = intersect_result
            end
        end
    end
    return scaffolds
end

function intersect(scaffold::Scaffold, bed_data::BedData)
    if !haskey(bed_data.scaffolds, scaffold.name)
        return nothing
    end

    return intersect(scaffold.features, bed_data.scaffolds[scaffold.name])
end

function intersect(genome::Genome, bed_data::BedData)
    scaffolds = Dict{String, IntervalMeta64}()
    for (scaffold_name, scaffold) in genome.scaffolds
        if haskey(bed_data.scaffolds, scaffold_name)
            intersect_result = intersect(scaffold, bed_data)
            if !isnothing(intersect_result)
                scaffolds[scaffold_name] = intersect_result
            end
        end
    end
    return scaffolds
end

"""
Merge the intervals of an interval tree into a sorted vector of disjoint,
closed `(start, end)` segments (1-based). Overlapping intervals are combined so
that each base is covered by at most one resulting segment.
"""
function merge_segments(tree::IntervalMeta64)
    segments = Tuple{Int, Int}[]
    ivs = sort!([(Int(iv.first), Int(iv.last)) for iv in tree])
    for (s, e) in ivs
        if !isempty(segments) && s <= segments[end][2]
            last_s, last_e = segments[end]
            segments[end] = (last_s, max(last_e, e))
        else
            push!(segments, (s, e))
        end
    end
    return segments
end

"""
Given a set of `BedData` measurements, compute how many of them cover each base
of the genome. The genome is returned as a dictionary keyed by scaffold name,
whose values are sparse arrays holding the per-base frequency (length = the
largest interval end seen on that scaffold).

When `merge` is `true` (the default), overlapping intervals *within a single
measurement* are merged first (via [`merge_segments`](@ref)), so each measurement
contributes at most 1 to a given base and the maximum possible value is the
number of measurements. Set `merge = false` to skip this step when the intervals
are already disjoint (e.g. ChIP-seq peak calls), in which case any within-measurement
overlaps will stack.

The element type is chosen to fit the measurement count: `UInt8` for up to 255
measurements and `UInt16` for up to 65535. More measurements raise an error.
"""
function calculate_frequency(measurements::Vector{BedData}; merge::Bool = true)
    n = length(measurements)
    T = if n <= typemax(UInt8)
        UInt8
    elseif n <= typemax(UInt16)
        UInt16
    else
        error("calculate_frequency supports at most $(Int(typemax(UInt16))) BedData measurements (received $n)")
    end

    # Every scaffold that appears in at least one measurement.
    scaffold_names = Set{String}()
    for measurement in measurements
        union!(scaffold_names, keys(measurement.scaffolds))
    end

    genome = Dict{String, SparseVector{T, Int}}()

    for name in scaffold_names
        # Difference array: +1 where a segment starts, -1 just past its end. The
        # running total while sweeping left-to-right is the per-base frequency
        # across measurements.
        deltas = Dict{Int, Int}()
        scaffold_len = 0

        for measurement in measurements
            haskey(measurement.scaffolds, name) || continue
            tree = measurement.scaffolds[name]
            segments = merge ? merge_segments(tree) :
                [(Int(iv.first), Int(iv.last)) for iv in tree]
            for (s, e) in segments
                deltas[s]     = get(deltas, s, 0) + 1
                deltas[e + 1] = get(deltas, e + 1, 0) - 1
                scaffold_len  = max(scaffold_len, e)
            end
        end

        isempty(deltas) && continue

        # Sweep the breakpoints in order, emitting a value for every covered base.
        breakpoints = sort!(collect(keys(deltas)))
        indices = Int[]
        values  = T[]
        coverage = 0
        for k in eachindex(breakpoints)
            p = breakpoints[k]
            coverage += deltas[p]
            if coverage > 0 && k < length(breakpoints)
                for base in p:(breakpoints[k + 1] - 1)
                    push!(indices, base)
                    push!(values, T(coverage))
                end
            end
        end

        genome[name] = sparsevec(indices, values, scaffold_len)
    end

    return genome
end

#= Base.show overloads =#

function Base.show(io::IO, b::BedData)
    n = length(b.scaffolds)
    print(io, "BedData($(n) scaffold$(n == 1 ? "" : "s"))")
end

function Base.show(io::IO, t::TabularData)
    r, c = size(t.table)
    print(io, "TabularData($(r)×$(c) $(eltype(t.table)))")
end

function Base.show(io::IO, b::BioSample)
    print(io, "BioSample(\"$(b.sample_id)\", tissue=$(b.tissue_type), species=\"$(b.species.name)\")")
end

Base.show(io::IO, a::AssayMethod) = print(io, "AssayMethod(\"$(a.name)\")")

function Base.show(io::IO, m::Measurement)
    print(io, "Measurement(\"$(basename(m.file_path))\", format=$(m.format))")
end

function Base.show(io::IO, a::Assay)
    print(io, "Assay(\"$(a.id)\", type=$(a.type))")
end

function Base.show(io::IO, s::Study)
    n = length(s.assays)
    print(io, "Study(\"$(s.id)\", \"$(s.title)\", $(s.date), $(n) assay$(n == 1 ? "" : "s"))")
end

export
    Study,
    Assay,
    AssayMethod,
    Measurement,
    BioSample,
    BedData,
    Data,
    TabularData,
    intersect,
    load_bed,
    load_table,
    merge_segments,
    calculate_frequency

end