module Data

using BED
using CSV
using CodecZlib
using DataFrames
using IntervalTrees
using ..BitCodes
using ..Reference

using Base: ImmutableDict

# Methylation is a sample-level signal like the rest of this module, but a
# positional one that needs no `Genome`, so it stays a self-contained submodule
# rather than being folded in here. Re-exported below, so both
# `BioinfoTools2.Data.Methylation` and `BioinfoTools2.Methylation` reach it.
include("data/methylation.jl")
using .Methylation

# Field layout of the 64-bit BED metadata code (see `pack_bed_code`).
const BED_STRAND_SHIFT = 32
const BED_STRAND_WIDTH = 8

"""
Interval tree type for BED data. Consists of a single sample's measurements,
associated with the `genome` the intervals are defined against.
"""
struct BedData
    genome::Genome
    scaffolds::Dict{String,IntervalTreeM64}
end

"""
Tabular data type. Consists of a matrix of numeric values from an arbitrary
number of samples, associated with the `genome` its samples are matched against.
"""
struct TabularData{T<:Number}
    genome::Genome
    variables::Vector{AbstractString}
    samples::Vector{Union{Nothing,Tuple{String,IntervalValue{UInt32,UInt64}}}}
    table::Matrix{T}
end

"""
A variable in an experiment.
Either it points to a set of covariates (child nodes) or contains data (leaf node).
The `value` of the variable is really the name of the incoming edge.
"""
struct Variable{T<:Number}
    value::AbstractString
    covariates::Union{Nothing,ImmutableDict{String,Variable{T}}}
    data::Union{Nothing,TabularData{T},BedData}

    function Variable{T}(
        value::AbstractString,
        covariates::Union{Nothing,ImmutableDict{String,Variable{T}}},
        data::Union{Nothing,TabularData{T},BedData},
    ) where {T<:Number}
        isnothing(covariates) ⊻ isnothing(data) ||
            throw(ArgumentError("Exactly one of `covariates` or `data` must be provided."))
        new{T}(value, covariates, data)
    end
end

# Recursively assemble the variable trie for an `Experiment`.
#
# `entries` pairs each sample's ordered variable-column values with its loaded
# data. `depth` is the 1-based variable column currently being grouped and
# `nvars` the total number of variable columns. Samples are grouped by their
# value at `depth` (first-seen order preserved); at the final column each group
# is a single validated sample and becomes a leaf `Variable` holding data,
# otherwise the group is recursed into to build that node's covariates.
function _build_variables(
    ::Type{T},
    entries::Vector{Tuple{Vector{String},D}},
    depth::Int,
    nvars::Int,
) where {T<:Number,D}
    order = String[]
    groups = Dict{String,Vector{Tuple{Vector{String},D}}}()
    for entry in entries
        key = entry[1][depth]
        haskey(groups, key) || push!(order, key)
        push!(get!(() -> Tuple{Vector{String},D}[], groups, key), entry)
    end

    variables = ImmutableDict{String,Variable{T}}()
    for key in order
        group = groups[key]
        node = if depth == nvars
            # A validated sample sheet guarantees one sample per leaf path.
            Variable{T}(key, nothing, group[1][2])
        else
            Variable{T}(key, _build_variables(T, group, depth + 1, nvars), nothing)
        end
        variables = ImmutableDict(variables, key => node)
    end

    return variables
end

# Load a single sample's data file, choosing the representation from the file
# extension: anything ending in `.bed` becomes `BedData`, everything else is
# read as a table and matched against `genome` as `TabularData{T}` (its numeric
# columns coerced to `T`).
function _load_sample(::Type{T}, genome::Genome, path::AbstractString) where {T<:Number}
    endswith(lowercase(path), ".bed") && return load_bed(genome, String(path))

    data_frame = CSV.read(path, DataFrame)
    for col in names(data_frame)[2:end]
        data_frame[!, col] = convert.(T, data_frame[!, col])
    end
    table = load_table(genome, data_frame)
    isnothing(table) && throw(ArgumentError("Could not load tabular data from \"$path\"."))
    return table
end

"""
Tree of experimental `Variable`s built from a `sample_sheet`, all associated
with `genome`.

The last column of `sample_sheet` must hold a path to each sample's data file;
every preceding column is treated as a categorical variable. Rows are organized
into a trie keyed by successive variable values, with each leaf holding the data
loaded from that sample's file. The data representation is chosen by file
extension: paths ending in `.bed` load as `BedData`, everything else loads as
`TabularData{T}` (numeric columns coerced to `T`). `T` defaults to `Float64`.
"""
struct Experiment{T<:Number}
    variables::ImmutableDict{String,Variable{T}}

    function Experiment{T}(genome::Genome, sample_sheet::DataFrame) where {T<:Number}
        ncols = ncol(sample_sheet)
        ncols >= 2 || throw(
            ArgumentError(
                "Sample sheet needs at least one variable column plus a file column.",
            ),
        )
        nvars = ncols - 1

        # 1. Every data file must be unique, and the variable columns must
        #    identify each sample uniquely (no duplicate variable-value paths).
        data_paths = sample_sheet[!, end]
        allunique(data_paths) ||
            throw(ArgumentError("File paths in the last column must be unique."))
        allunique(Tuple(row) for row in eachrow(sample_sheet[!, 1:nvars])) || throw(
            ArgumentError(
                "Each row must have a unique combination of variable-column values.",
            ),
        )

        # 2. Every purported file must exist.
        all(isfile, data_paths) || throw(
            ArgumentError(
                "All paths in the last column of the sample sheet must be valid file paths.",
            ),
        )

        # 3. Load each sample's data concurrently (representation inferred from
        #    the file extension), pairing it with that row's variable values
        #    (kept as strings regardless of the column type).
        D = Union{TabularData{T},BedData}
        tasks = map(eachrow(sample_sheet)) do row
            path = String(row[end])
            variable_values = String[string(row[c]) for c = 1:nvars]
            Threads.@spawn (variable_values, _load_sample(T, genome, path))
        end
        entries = Tuple{Vector{String},D}[fetch(task) for task in tasks]

        # 4. Fold the loaded samples into the immutable variable trie.
        return new{T}(_build_variables(T, entries, 1, nvars))
    end
end

Experiment(genome::Genome, sample_sheet::DataFrame) =
    Experiment{Float64}(genome, sample_sheet)

#= Methods =#

"""
Match a data frame's rows against `genome` by ID, returning a `TabularData`.

The first column must hold sample IDs (e.g. gene or transcript accessions) and
every other column must be numeric; `nothing` is returned (with a warning) when
that does not hold. Rows whose ID matches no feature keep a `nothing` sample
slot. The whole table is resolved in a single O(features) genome walk, which
stops early once every row has matched.
"""
function load_table(genome::Genome, data_frame::DataFrame)
    sample_col_correct = eltype(data_frame[!, 1]) <: AbstractString
    data_cols_correct = all(eltype(data_frame[!, i]) <: Number for i = 2:ncol(data_frame))
    sample_col_correct || @warn "Sample column incorrect format"
    data_cols_correct || @warn "Data columns incorrect format"
    (sample_col_correct && data_cols_correct) || return nothing

    n_rows = nrow(data_frame)

    # Output vector, parallel to the table rows. A slot stays `nothing` if its
    # sample ID never matches any feature's metadata ID. A match stores the
    # scaffold name and the full interval (start/end coords + 64-bit metadata
    # code), so the feature type can be filtered out of these later without
    # re-reading the table.
    samples =
        Vector{Union{Nothing,Tuple{String,IntervalValue{UInt32,UInt64}}}}(nothing, n_rows)

    # Samples still awaiting a match: sample_id => original row index. Popping a
    # match makes the whole search O(features) with an early exit.
    remaining = Dict{String,Int}()
    sizehint!(remaining, n_rows)
    for (row, id) in enumerate(data_frame[!, 1])
        remaining[String(id)] = row
    end

    for (scaffold_name, scaffold) in genome.scaffolds
        isempty(remaining) && break

        for interval in scaffold.features
            id = Reference.get_metadata_id(genome, Reference.parse_index(interval.value))
            isnothing(id) && continue

            row = get(remaining, id, nothing)
            if !isnothing(row)
                samples[row] = (scaffold_name, interval)
                delete!(remaining, id)
                isempty(remaining) && break
            end
        end
    end

    return TabularData(
        genome,
        Vector{AbstractString}(names(data_frame)[2:end]),
        samples,
        Matrix(data_frame[!, 2:end]),
    )
end

#= TabularData indexing =#

# Integer / cartesian indexing behaves exactly like indexing the underlying
# matrix (e.g. `tab[2, 2]`, `tab[1:10, :]`).
Base.getindex(t::TabularData, inds...) = getindex(t.table, inds...)

"""
Resolve the ID of a single sample by looking its metadata up in `genome`.
Returns `nothing` for empty slots or features without metadata.
"""
function sample_id(
    genome::Genome,
    sample::Union{Nothing,Tuple{String,IntervalValue{UInt32,UInt64}}},
)
    isnothing(sample) && return nothing
    return Reference.get_metadata_id(genome, Reference.parse_index(sample[2].value))
end

# Rows of `t` whose sample ID satisfies `matches`, in table order.
_matching_rows(t::TabularData, matches::F) where {F} =
    findall(sample -> matches(sample_id(t.genome, sample)), t.samples)

# `t` restricted to `rows`, sharing nothing mutable with the original.
_subset(t::TabularData, rows::AbstractVector{Int}) =
    TabularData(t.genome, copy(t.variables), t.samples[rows], t.table[rows, :])

"""
Look up a single sample by its (putative) ID string, returning the first match
as a one-row `TabularData`, or `nothing` when no sample matches.

This is intentionally O(samples) per lookup (it resolves each sample's ID on
demand) rather than maintaining an ID index.
"""
function Base.getindex(t::TabularData, id::AbstractString)
    row = findfirst(sample -> sample_id(t.genome, sample) == id, t.samples)
    return isnothing(row) ? nothing : _subset(t, [row])
end

"""
Look up every sample whose (putative) ID is contained in `ids`, returning a
`TabularData` sub-table (with matching `variables`, `samples` and `table`
rows). IDs with no matching sample are silently skipped; row order follows the
original table.
"""
function Base.getindex(t::TabularData, ids::AbstractVector{<:AbstractString})
    wanted = Set{String}(ids)
    return _subset(t, _matching_rows(t, id -> !isnothing(id) && id in wanted))
end

"""
Select a subset of a `TabularData`'s variables (columns) by name, returning a
new `TabularData` restricted to the matching variables (or `nothing` when none
of the requested variables are present).
"""
function DataFrames.select(
    t::TabularData,
    variables::Union{AbstractString,AbstractVector{<:AbstractString}},
)::Union{Nothing,TabularData}
    wanted = variables isa AbstractString ? [variables] : variables
    columns = findall(in(wanted), t.variables)
    isempty(columns) && return nothing
    return TabularData(t.genome, t.variables[columns], t.samples, t.table[:, columns])
end

"""
Convert a `TabularData` back into a `DataFrame`. The first column, `ID`, holds
the ID that each sample's metadata index points to (resolved from the table's
own genome; unmatched samples become `missing`), followed by one column per
entry in `variables` carrying the corresponding `table` column.
"""
function DataFrames.DataFrame(t::TabularData)
    data_frame = DataFrame()
    data_frame[!, :ID] =
        Union{Missing,String}[something(sample_id(t.genome, s), missing) for s in t.samples]
    for (column, variable) in enumerate(t.variables)
        data_frame[!, string(variable)] = t.table[:, column]
    end
    return data_frame
end

#= BED data =#

"""
Pack a BED interval's strand into a 64-bit metadata code. Bits 33-40 hold the
package-wide 2-bit strand code (see `BitCodes`); every other bit is reserved.

|  64-41  | 40-33  |   32-1   |
|---------|--------|----------|
| <NULL>  | Strand |  <NULL>  |
"""
pack_bed_code(strand::UInt8) =
    set_field(UInt64(0), strand, BED_STRAND_SHIFT, BED_STRAND_WIDTH)

"""Extract the strand code from a BED interval's metadata code."""
parse_bed_strand(code::UInt64) = UInt8(get_field(code, BED_STRAND_SHIFT, BED_STRAND_WIDTH))

"""
Convert a BED record's strand field to the package-wide 2-bit strand code
(see `BitCodes`). Records without a strand field are encoded as `STRAND_NA`.
"""
function bed_record_strand(record::BED.Record)
    BED.hasstrand(record) || return STRAND_NA
    return strand_code(BED.strand(record))
end

"""
Load a BED file (optionally gzipped) into a [`BedData`](@ref) associated with
`genome`. BED's 0-based half-open coordinates are converted to the 1-based
closed intervals used by gene features.
"""
function load_bed(genome::Genome, file_path::String)
    scaffolds = Dict{String,IntervalTreeM64}()

    open(file_path) do fh
        reader =
            endswith(file_path, ".gz") ? BED.Reader(GzipDecompressorStream(fh)) :
            BED.Reader(fh)

        record = BED.Record()
        try
            while !eof(reader)
                empty!(record)
                read!(reader, record)
                # Process if filled: handles both normal reads and the last record
                # in files without a trailing newline (EOFError thrown after fill).
                BED.isfilled(record) || continue

                code = pack_bed_code(bed_record_strand(record))
                tree = get!(IntervalTreeM64, scaffolds, BED.chrom(record))
                push!(
                    tree,
                    IntervalValue(
                        UInt32(BED.chromstart(record)),
                        UInt32(BED.chromend(record)),
                        code,
                    ),
                )
            end
        finally
            close(reader)
        end
    end

    return BedData(genome, scaffolds)
end

#= Set operations on interval trees =#

"""
Given two `IntervalTreeM64` trees 'a' and 'b', return the tree representing their intersection,
with the metadata from tree 'a' kept as the metadata for the intersection.

**NOTE:** multiple intervals from tree 'b' can intersect one interval from tree 'a', and
therefore multiple intervals in the return tree can have the same metadata (the same feature can
be present in the return tree multiple times in fragments).
"""
function intersect(tree_a::IntervalTreeM64, tree_b::IntervalTreeM64)::IntervalTreeM64
    intersection = IntervalTreeM64()
    for (left, right) in IntervalTrees.intersect(tree_a, tree_b)
        push!(
            intersection,
            IntervalValue(
                max(left.first, right.first),
                min(left.last, right.last),
                left.value,
            ),
        )
    end
    return intersection
end

"""
Intersect the `feature`-type features of `scaffold` with the matching scaffold
of `bed_data`. Returns `nothing` when the scaffold is absent from `bed_data` or
`feature` is not a known SO term.
"""
function intersect(
    scaffold::Scaffold,
    bed_data::BedData,
    feature::Union{AbstractString,Symbol},
)
    tree = get(bed_data.scaffolds, scaffold.name, nothing)
    isnothing(tree) && return nothing
    feature_intervals = get_feature(scaffold, feature)
    isnothing(feature_intervals) && return nothing
    return intersect(feature_intervals, tree)
end

"""
Intersect every feature of `scaffold` with the matching scaffold of `bed_data`.
Returns `nothing` when the scaffold is absent from `bed_data`.
"""
function intersect(scaffold::Scaffold, bed_data::BedData)::Union{Nothing,IntervalTreeM64}
    tree = get(bed_data.scaffolds, scaffold.name, nothing)
    isnothing(tree) && return nothing
    return intersect(scaffold.features, tree)
end

"""
Intersect `genome` with `bed_data` scaffold by scaffold, restricted to
`feature`-type features when one is given. Scaffolds with no intersection are
omitted from the result.
"""
function intersect(genome::Genome, bed_data::BedData, args...)::Dict{String,IntervalTreeM64}
    scaffolds = Dict{String,IntervalTreeM64}()
    for (scaffold_name, scaffold) in genome.scaffolds
        result = intersect(scaffold, bed_data, args...)
        isnothing(result) || (scaffolds[scaffold_name] = result)
    end
    return scaffolds
end

"""
Left-join two interval trees, preserving every interval of `treeL` and pairing
each with a (possibly `nothing`) interval from `treeR`.

The returned value is a lazy iterator of `(left, right)` tuples where `left` is
always an interval from `treeL` and `right` is either a matching interval from
`treeR` or `nothing`. Left intervals are visited in the tree's natural (sorted)
order; a left interval with several right matches yields one tuple per match,
while a left interval with no match yields a single `(left, nothing)` tuple.
Right intervals that never match a left interval are dropped.

How intervals are matched is controlled by `on`:

- `:metadata` (default) : match on the 64-bit metadata `value`.
- `:start`              : match on the start position (`first`).
- `:end`                : match on the end position (`last`).
- `:interval`           : match only when both `first` and `last` are equal.
"""
function leftjoin(treeL::IntervalTreeM64, treeR::IntervalTreeM64, on::Symbol = :metadata)
    on === :metadata && return _leftjoin(treeL, treeR, iv -> iv.value, UInt64)
    on === :start && return _leftjoin(treeL, treeR, iv -> iv.first, UInt32)
    on === :end && return _leftjoin(treeL, treeR, iv -> iv.last, UInt32)
    on === :interval &&
        return _leftjoin(treeL, treeR, iv -> (iv.first, iv.last), Tuple{UInt32,UInt32})
    throw(
        ArgumentError("`on` must be one of :metadata, :start, :end, :interval (got :$on)"),
    )
end

# Left-join `treeL` against `treeR` on `keyfn`, whose return type is `KT`.
function _leftjoin(
    treeL::IntervalTreeM64,
    treeR::IntervalTreeM64,
    keyfn::F,
    ::Type{KT},
) where {F,KT}
    IV = eltype(treeL)   # IntervalValue{UInt32, UInt64}

    # Index the right tree by join key: key => right intervals sharing that key.
    right_index = Dict{KT,Vector{IV}}()
    for right in treeR
        push!(get!(() -> IV[], right_index, keyfn(right)), right)
    end

    # For each left interval, emit one tuple per matching right interval, or a
    # single `(left, nothing)` tuple when there is no match. `flatten` keeps the
    # whole thing lazy.
    return Iterators.flatten(
        let matches = get(right_index, keyfn(left), nothing)
            isnothing(matches) ? ((left, nothing),) : ((left, right) for right in matches)
        end for left in treeL
    )
end

"""
Merge the intervals of an interval tree into a sorted vector of disjoint,
closed `(start, end)` segments (1-based). Overlapping intervals are combined so
that each base is covered by at most one resulting segment.
"""
function merge_segments(tree::IntervalTreeM64)
    segments = Tuple{Int,Int}[]
    for (start_pos, end_pos) in sort!([(Int(iv.first), Int(iv.last)) for iv in tree])
        if !isempty(segments) && start_pos <= segments[end][2]
            segments[end] = (segments[end][1], max(segments[end][2], end_pos))
        else
            push!(segments, (start_pos, end_pos))
        end
    end
    return segments
end

#= Base.show overloads =#

function Base.show(io::IO, b::BedData)
    n = length(b.scaffolds)
    print(io, "BedData($(n) scaffold$(n == 1 ? "" : "s"))")
end

function Base.show(io::IO, t::TabularData)
    rows, columns = size(t.table)
    print(io, "TabularData($(rows)×$(columns) $(eltype(t.table)))")
end

export BedData,
    Methylation, TabularData, intersect, leftjoin, load_bed, load_table, merge_segments

end
