module Exploration

# Scoped, not a blanket `using DataFrames` — that would also pull in
# `DataFrames.leftjoin`, colliding with `Data.leftjoin` used by `coverage`.
using DataFrames: DataFrame, nrow
using Interpolations
using KernelDensity
using SparseArrays
using StatsBase

using ..Reference
using ..Data

# Scoped for the same reason as the DataFrames import above: `Methylation`
# exports a large surface (`merge_calls`, `load_bismark`, …) that this module
# has no use for.
using ..Data.Methylation:
    MethylationData, CTX_CPG, find_calls_in_range, get_context, get_depth, meth_fraction

#= Feature regions =#

"""
Given one feature's interval, work out the region a profile covers. Returns
`(feature_id, negative, region_start, region_end)`, or `nothing` when the
feature's metadata ID does not resolve.

`clip` holds `region_start` at 1: the methylation path leaves it unclipped, so
that index `flank + 1` is the feature's first base for every feature.
"""
function _feature_region(genome::Genome, interval, flank::Integer; clip::Bool)
    found = Reference.feature_id(genome, interval)
    isnothing(found) && return nothing

    region_start = Int(interval.first) - Int(flank)
    clip && (region_start = max(1, region_start))
    return (
        found,
        Reference.parse_strand(interval.value) == get_strand('-'),
        region_start,
        Int(interval.last) + Int(flank),
    )
end

"""
Given a genomic base, map it into the region. Returns its 1-based index there,
counted from the region's far end for a `negative`-strand feature so that index
1 is always the feature's 5' end.
"""
@inline _region_index(base::Integer, region_start::Int, region_end::Int, negative::Bool) =
    negative ? region_end - Int(base) + 1 : Int(base) - region_start + 1

#= Per-feature coverage =#

"""
    coverage(data::BedData, feature; filter_zeros = false)

Given BED intervals, score how much of each `feature`-type feature of
`data.genome` they cover. Returns a `Dict` mapping scaffold name to its vector of
`covered_length / feature_length` fractions in `[0, 1]`, in the scaffold's own
feature order; scaffolds absent from `data` are omitted.

Set `filter_zeros = true` to drop the uncovered (`0.0`) entries.
"""
function coverage(
    data::BedData,
    feature::Union{AbstractString,Symbol};
    filter_zeros::Bool = false,
)::Dict{String,Vector{Float64}}
    intersection = Data.intersect(data.genome, data, feature)
    scaffolds = Dict{String,Vector{Float64}}()

    for (scaffold_name, scaffold) in data.genome.scaffolds
        haskey(intersection, scaffold_name) || continue

        feature_tree = scaffold.features
        fractions = Float64[]
        sizehint!(fractions, length(feature_tree))
        for (whole, covered) in leftjoin(feature_tree, intersection[scaffold_name])
            if isnothing(covered)
                push!(fractions, 0.0)
            else
                feature_length = whole.last - whole.first + 1
                covered_length = covered.last - covered.first + 1
                @assert feature_length >= covered_length
                push!(fractions, covered_length / feature_length)
            end
        end

        filter_zeros && filter!(!iszero, fractions)
        scaffolds[scaffold_name] = fractions
    end
    return scaffolds
end

"""
    coverage(genome, data::MethylationData, feature; filter_zeros = false)

Given single-base methylation calls, score each `feature`-type feature of
`genome` by them. Returns a `Dict` mapping scaffold name to its vector of mean
per-base levels in `[0, 1]`, on the same scale as the [`BedData`](@ref) method;
scaffolds absent from `data` are omitted.

**NOTE:** a feature's score is the summed methylation fraction over its whole
length, so **bases with no call count as zero** and the score reflects
methylation level and cytosine density together. Strand and context are ignored.
"""
function coverage(
    genome::Genome,
    data::MethylationData,
    feature::Union{AbstractString,Symbol};
    filter_zeros::Bool = false,
)::Dict{String,Vector{Float64}}
    scaffolds = Dict{String,Vector{Float64}}()

    for (scaffold_name, tree) in get_feature(genome, feature)
        haskey(data, scaffold_name) || continue
        calls = data[scaffold_name]

        mean_levels = Float64[]
        sizehint!(mean_levels, length(tree))
        for interval in tree
            total_fraction = 0.0
            for call in find_calls_in_range(calls, interval.first, interval.last)
                fraction = meth_fraction(call)
                # `meth_fraction` is NaN for a site with no coverage at all.
                isnan(fraction) || (total_fraction += fraction)
            end
            push!(
                mean_levels,
                total_fraction / (Int(interval.last) - Int(interval.first) + 1),
            )
        end

        filter_zeros && filter!(!iszero, mean_levels)
        scaffolds[scaffold_name] = mean_levels
    end
    return scaffolds
end

#= Density estimates =#

"""
Given per-scaffold score vectors, fit a kernel density estimate to each.
Returns a `Dict` mapping scaffold name to its `UnivariateKDE`, or to `nothing`
where there was nothing to fit.
"""
_kde_by_scaffold(by_scaffold::Dict{String,Vector{Float64}}) =
    Dict{String,Union{Nothing,UnivariateKDE}}(
        name => isempty(scores) ? nothing : KernelDensity.kde(scores) for
        (name, scores) in by_scaffold
    )

"""
    kde(data::BedData, feature; filter_zeros = false)

Given BED intervals, estimate the distribution of the per-feature coverage
fractions [`coverage`](@ref) scores them at. Returns a `Dict` mapping scaffold
name to a `UnivariateKDE`, or to `nothing` where there was nothing to fit.
"""
kde(data::BedData, feature::Union{AbstractString,Symbol}; filter_zeros::Bool = false) =
    _kde_by_scaffold(coverage(data, feature; filter_zeros = filter_zeros))

"""
    kde(genome, data::MethylationData, feature; filter_zeros = false)

Given methylation calls, estimate the distribution of the per-feature mean
levels [`coverage`](@ref) scores them at. Returns a `Dict` mapping scaffold name
to a `UnivariateKDE`, or to `nothing` where there was nothing to fit.
"""
kde(
    genome::Genome,
    data::MethylationData,
    feature::Union{AbstractString,Symbol};
    filter_zeros::Bool = false,
) = _kde_by_scaffold(coverage(genome, data, feature; filter_zeros = filter_zeros))

"""
    kde(data::TabularData; filter_zeros = false, transform = identity)

Given a table, estimate the distribution of every value in it after applying
`transform`. Returns a `UnivariateKDE`, or `nothing` when nothing is left to
fit. Non-finite values are dropped, as are zeros when `filter_zeros = true`.
"""
function kde(data::TabularData; filter_zeros::Bool = false, transform::Function = identity)
    flat = [transform(value) for value in data.table]
    filter!(value -> isfinite(value) && !(filter_zeros && iszero(value)), flat)
    return isempty(flat) ? nothing : KernelDensity.kde(flat)
end

#= Quantile binning =#

"""
Given the sorted quantile `edges` (length `n_bins + 1`), bin `value`. Returns
its 1-based bin; a value on an edge falls into the lower bin, and anything at or
above the top edge into the last one.
"""
_quantile_bin(edges::AbstractVector, value::Real, n_bins::Int) =
    min(searchsortedfirst(view(edges, 2:lastindex(edges)), value), n_bins)

"""
Given the arguments the rank-based `quantiles` methods share, validate them.
Returns `nothing`; `noun` names what `ranking` selects, for the error message.
"""
function _check_ranking(n_bins::Int, ranking::Vector{String}, noun::String)
    n_bins >= 1 ||
        throw(ArgumentError("`quantiles` must be a positive integer (got $n_bins)"))
    isempty(ranking) && throw(ArgumentError("`ranking` must name at least one $noun"))
    return nothing
end

"""
Given ranking keys — tuples of ranking values with the original row index last —
sort them and cut the ranking into `n_bins` equal-sized (up to an off-by-one)
bins. Returns `(row, bin)` pairs in ranked order.
"""
function _ranked_bins(ranking_keys::Vector, n_bins::Int)
    sort!(ranking_keys)
    n_rows = length(ranking_keys)
    return [
        (key[end]::Int, cld(rank_pos * n_bins, n_rows)) for
        (rank_pos, key) in enumerate(ranking_keys)
    ]
end

"""
    quantiles(data::TabularData; quantiles = 4, merge = mean)

Given a table, bin each sample by a scalar summary of its row: `merge` (default
`mean`) collapses the row, and the merged values define `quantiles + 1` quantile
edges. Returns `(FeatureRecord, merged_value, quantile_index)` tuples in the
sample order of `data`, bin `1` holding the lowest values.

Each `FeatureRecord` is looked up in `data.genome` by the sample's 32-bit
metadata index; unmatched samples and unresolvable features are skipped.
"""
function quantiles(data::TabularData; quantiles::Int = 4, merge = mean)
    quantiles >= 1 ||
        throw(ArgumentError("`quantiles` must be a positive integer (got $quantiles)"))

    # Collapse each matched sample's row to a scalar, remembering its metadata index.
    indices = UInt32[]
    merged = Float64[]
    for (row, sample) in enumerate(data.samples)
        isnothing(sample) && continue
        push!(indices, Reference.parse_index(sample[2].value))
        push!(merged, Float64(merge(data.table[row, :])))
    end

    result = Tuple{FeatureRecord,Float64,Int}[]
    isempty(merged) && return result

    # `quantiles + 1` edges spanning the observed range of merged values.
    edges = quantile(merged, range(0, 1; length = quantiles + 1))

    # One genome walk for every index at once; a lookup per index would be
    # O(rows * features).
    records = data.genome[indices]

    for (meta_index, value) in zip(indices, merged)
        record = get(records, meta_index, nothing)
        isnothing(record) && continue
        push!(result, (record, value, _quantile_bin(edges, value, quantiles)))
    end

    return result
end

"""
    quantiles(data::TabularData, ranking::Vector{String}; quantiles::Int = 4)

Given a table and the variables to rank its samples by — the first primary, each
later one breaking ties in the one before it, original sample order last — cut
that ranking into `quantiles` equal-sized (up to an off-by-one) bins. Returns
`(FeatureRecord, quantile_index)` tuples in ranked order.

**NOTE:** samples whose feature does not resolve are skipped from the result but
still occupy a rank position, and so still affect bin sizing.
"""
function quantiles(data::TabularData, ranking::Vector{String}; quantiles::Int = 4)
    _check_ranking(quantiles, ranking, "variable")

    columns = map(ranking) do name
        column = findfirst(==(name), data.variables)
        isnothing(column) && throw(ArgumentError("Unknown ranking variable: \"$name\""))
        column
    end

    # Matched samples only; the row index itself is the final tie-break.
    rows = findall(!isnothing, data.samples)
    isempty(rows) && return Tuple{FeatureRecord,Int}[]
    ranked = _ranked_bins(
        [(Tuple(data.table[row, c] for c in columns)..., row) for row in rows],
        quantiles,
    )

    # Resolve every ranked sample's feature in one genome walk; a lookup per rank
    # position would be O(rows * features).
    meta_indices =
        UInt32[Reference.parse_index(data.samples[row][2].value) for (row, _) in ranked]
    records = data.genome[meta_indices]

    result = Tuple{FeatureRecord,Int}[]
    sizehint!(result, length(ranked))
    for (meta_index, (_, bin)) in zip(meta_indices, ranked)
        record = get(records, meta_index, nothing)
        isnothing(record) || push!(result, (record, bin))
    end
    return result
end

"""
    quantiles(pairs::DataFrame, ranking::Vector{String}; quantiles::Int = 4) -> DataFrame

Given a paralog-pair table (shaped like `ParalogGroup`'s constructor input, or
`rbh`'s output) and the columns to rank its rows by — the first primary, each
later one breaking ties in the one before it, original row order last — cut that
ranking into `quantiles` equal-sized (up to an off-by-one) bins. Returns `pairs`
with a 1-based `"quantile"` column appended; row order is unchanged.
"""
function quantiles(pairs::DataFrame, ranking::Vector{String}; quantiles::Int = 4)
    _check_ranking(quantiles, ranking, "column")
    for name in ranking
        name in names(pairs) || throw(ArgumentError("Unknown ranking column: \"$name\""))
    end

    n_rows = nrow(pairs)
    bins = Vector{Int}(undef, n_rows)
    if n_rows > 0
        ranked = _ranked_bins(
            [(Tuple(pairs[row, name] for name in ranking)..., row) for row = 1:n_rows],
            quantiles,
        )
        for (row, bin) in ranked
            bins[row] = bin
        end
    end

    result = copy(pairs)
    result.quantile = bins
    return result
end

#= Genome-wide coverage frequency =#

"""
Given a measurement count, pick the narrowest element type that can hold it.
Returns `UInt8` up to 255 measurements and `UInt16` up to 65535; more is an
error.
"""
function _frequency_eltype(n_measurements::Int)
    n_measurements <= typemax(UInt8) && return UInt8
    n_measurements <= typemax(UInt16) && return UInt16
    error(
        "calculate_frequency supports at most $(Int(typemax(UInt16))) BedData measurements (received $n_measurements)",
    )
end

"""
Given a set of measurements, list every scaffold at least one of them covers.
Returns the names as a `Vector`, which `Threads.@threads` needs in place of a
`Set`.
"""
function _measured_scaffolds(measurements::Vector{BedData})
    scaffold_names = String[]
    seen = Set{String}()
    for measurement in measurements, name in keys(measurement.scaffolds)
        name in seen || (push!(seen, name); push!(scaffold_names, name))
    end
    return scaffold_names
end

"""
Given one scaffold, collect the measurements' segment boundaries over it.
Returns `(events, scaffold_length)`, the events a sorted difference array of
`(position, delta)`: `+1` where a segment starts, `-1` just past its end.

`merge` merges the intervals within a measurement first, so that it contributes
at most 1 to a base.
"""
function _coverage_events(measurements::Vector{BedData}, name::String, merge::Bool)
    events = Tuple{Int,Int}[]
    scaffold_length = 0

    for measurement in measurements
        tree = get(measurement.scaffolds, name, nothing)
        isnothing(tree) && continue
        segments =
            merge ? merge_segments(tree) : [(Int(iv.first), Int(iv.last)) for iv in tree]
        for (start_pos, end_pos) in segments
            push!(events, (start_pos, 1))
            push!(events, (end_pos + 1, -1))
            scaffold_length = max(scaffold_length, end_pos)
        end
    end
    return sort!(events), scaffold_length
end

"""
Given a scaffold's difference array, sweep its breakpoints in order. Returns a
`SparseVector{T}` holding, for every covered base, how many measurements cover
it.
"""
function _sweep_events(
    events::Vector{Tuple{Int,Int}},
    scaffold_length::Int,
    ::Type{T},
) where {T}
    indices = Int[]
    counts = T[]
    depth = 0
    event = 1
    n_events = length(events)

    while event <= n_events
        position = events[event][1]
        while event <= n_events && events[event][1] == position
            depth += events[event][2]
            event += 1
        end
        # The final breakpoint always closes the last segment (depth 0).
        if depth > 0 && event <= n_events
            for base = position:(events[event][1]-1)
                push!(indices, base)
                push!(counts, T(depth))
            end
        end
    end
    return sparsevec(indices, counts, scaffold_length)
end

"""
    calculate_frequency(measurements; merge = true)

Given a set of `BedData` measurements, count how many of them cover each base of
the genome. Returns a `Dict` mapping scaffold name to a `SparseVector` of
per-base frequencies, as long as the largest interval end seen on that scaffold.

`merge = true` (the default) merges the overlapping intervals *within* one
measurement first, capping the result at the measurement count; set it to
`false` for already-disjoint intervals (e.g. ChIP-seq peak calls), which then
stack.
"""
function calculate_frequency(measurements::Vector{BedData}; merge::Bool = true)
    T = _frequency_eltype(length(measurements))
    scaffold_names = _measured_scaffolds(measurements)

    # Pre-populate every key so the parallel loop only overwrites existing
    # entries. Inserting new keys concurrently would race on the Dict's internal
    # structure; overwriting the value of an existing key does not.
    genome =
        Dict{String,SparseVector{T,Int}}(name => spzeros(T, 0) for name in scaffold_names)

    Threads.@threads for name in scaffold_names
        events, scaffold_length = _coverage_events(measurements, name, merge)
        isempty(events) || (genome[name] = _sweep_events(events, scaffold_length, T))
    end

    return genome
end

"""
Per-feature, per-base coverage counts across a set of `BedData` measurements,
together with the number of measurements (`n`) they were computed from.

`features` maps a feature ID to a `SparseVector{UInt32}` of raw overlap counts —
one entry per base of the (flanked) feature, oriented in the direction of
transcription (index 1 is the feature's 5' end). Dividing a count by `n` gives
the fraction of measurements covering that base.
"""
struct FeatureFrequency
    n::Int
    features::Dict{String,SparseVector{UInt32,Int}}
end

"""
    feature_frequency(genome, feature, frequency, n; flank = 500)

Given a per-base `frequency` dictionary (as [`calculate_frequency`](@ref)
returns), project it onto every `feature`-type feature of `genome`. Returns a
[`FeatureFrequency`](@ref) carrying the measurement count `n`, so that per-base
frequencies can be recovered by division.

Each feature holds the padded region `[first - flank, last + flank]`, re-indexed
from 1 and reversed on the negative strand so index 1 is its 5' end. Features
whose metadata ID cannot be resolved are skipped.
"""
function feature_frequency(
    genome::Genome,
    feature::Union{AbstractString,Symbol},
    frequency::AbstractDict{String,<:AbstractVector},
    n::Integer;
    flank::Integer = 500,
)
    features = Dict{String,SparseVector{UInt32,Int}}()

    for (scaffold_name, tree) in get_feature(genome, feature)
        counts = get(frequency, scaffold_name, nothing)
        # Nonzero (base, count) pairs for this scaffold, ascending by position.
        base_indices, base_counts = isnothing(counts) ? (Int[], UInt32[]) : findnz(counts)

        for interval in tree
            region = _feature_region(genome, interval, flank; clip = true)
            isnothing(region) && continue
            feature_id, negative, region_start, region_end = region

            # Slice of nonzero bases falling inside the padded region.
            entries =
                searchsortedfirst(base_indices, region_start):searchsortedlast(
                    base_indices,
                    region_end,
                )
            features[feature_id] = sparsevec(
                Int[
                    _region_index(base_indices[e], region_start, region_end, negative)
                    for e in entries
                ],
                UInt32[base_counts[e] for e in entries],
                region_end - region_start + 1,
            )
        end
    end

    return FeatureFrequency(n, features)
end

#= Per-feature methylation levels =#

"""
Minimum read depth a call must have to enter a methylation profile. Shallow
sites carry a level that is mostly sampling noise — a 1-read site reads as 0% or
100% — so they are dropped rather than weighted down.
"""
const DEFAULT_MIN_DEPTH = UInt32(5)

"""
    weight_of(depth, weight_transform)

Given a read depth, turn it into the statistical weight the profile functions
give it. Returns `weight_transform(depth)` as a `Float64`, throwing an
`ArgumentError` for a weight that is negative or not finite.

**NOTE:** a transform may legitimately return `0` — `log(1) == 0` — which drops
that observation; prefer `log1p` over `log` unless that is what you want.
"""
@inline function weight_of(depth::Real, weight_transform)
    weight = Float64(weight_transform(Float64(depth)))
    (isfinite(weight) && weight >= 0) || throw(
        ArgumentError(
            "`weight_transform` returned $weight for depth $depth; weights must be finite and non-negative",
        ),
    )
    return weight
end

"""
One feature's per-base methylation, oriented in the direction of transcription
(index 1 is the 5' end of the flanked region).

- `levels`: depth-weighted methylation fraction in `[0, 1]` at each base.
- `weights`: total read depth at each base. **A nonzero weight is what marks a
  base as measured**, since a base measured at 0% methylation and a base with no
  cytosine are both structural zeros in `levels`.
"""
struct FeatureLevels
    levels::SparseVector{Float32,Int}
    weights::SparseVector{UInt32,Int}
end

"""
Per-feature, per-base methylation levels, together with the filters they were
built under.

The methylation counterpart of [`FeatureFrequency`](@ref), kept separate because
the two carry different quantities: `FeatureFrequency` holds overlap *counts* to
be divided by a measurement count, whereas this holds a *level* per base plus
the depth it rests on, and must distinguish an unmeasured base from one measured
at zero.
"""
struct MethylationFrequency
    min_depth::UInt32
    context::Union{Nothing,UInt8}
    features::Dict{String,FeatureLevels}
end

"""
Given the calls over one region, pool them per base. Returns
`(indices, levels, weights)` in region coordinates, several calls on one base —
the two strands of a CpG, say — combining in proportion to their depths.

Calls shallower than `min_depth`, or outside `context` when one is given, are
dropped.
"""
function _region_levels(
    region_calls,
    region_start::Int,
    region_end::Int,
    negative::Bool,
    min_depth::UInt32,
    context::Union{Nothing,UInt8},
)
    indices = Int[]
    level_values = Float32[]
    weight_values = UInt32[]

    # Calls are position-sorted, so the (rare) several calls sharing a base form
    # one contiguous run.
    index = 1
    n_calls = length(region_calls)
    while index <= n_calls
        position = region_calls.pos[index]
        weighted_fraction = 0.0
        depth_total = 0

        while index <= n_calls && region_calls.pos[index] == position
            call = region_calls[index]
            index += 1
            depth = get_depth(call)
            depth < min_depth && continue
            isnothing(context) || get_context(call) == context || continue
            fraction = meth_fraction(call)
            isnan(fraction) && continue
            # Always linear in depth: this pools reads at one cytosine back into
            # that cytosine's own methylation fraction, which only comes out
            # right weighted by read count.
            weighted_fraction += fraction * Int(depth)
            depth_total += Int(depth)
        end
        depth_total == 0 && continue

        push!(indices, _region_index(position, region_start, region_end, negative))
        push!(level_values, Float32(weighted_fraction / depth_total))
        push!(weight_values, UInt32(min(depth_total, Int(typemax(UInt32)))))
    end
    return indices, level_values, weight_values
end

"""
    feature_frequency(genome, feature, data::MethylationData; flank = 500,
                      min_depth = DEFAULT_MIN_DEPTH, context = CTX_CPG)

Given single-base methylation calls, project them onto every `feature`-type
feature of `genome`. Returns a [`MethylationFrequency`](@ref) holding each
feature's padded region `[first - flank, last + flank]`, re-indexed from 1 and
reversed on the negative strand.

**NOTE:** unlike the [`BedData`](@ref) pipeline the region is *not* clipped at
the scaffold start, so index `flank + 1` is the feature's first base for every
feature.
"""
function feature_frequency(
    genome::Genome,
    feature::Union{AbstractString,Symbol},
    data::MethylationData;
    flank::Integer = 500,
    min_depth::Integer = DEFAULT_MIN_DEPTH,
    context::Union{Nothing,Integer} = CTX_CPG,
)
    min_depth_bits = UInt32(min_depth)
    context_bits = isnothing(context) ? nothing : UInt8(context)
    features = Dict{String,FeatureLevels}()

    for (scaffold_name, tree) in get_feature(genome, feature)
        haskey(data, scaffold_name) || continue
        calls = data[scaffold_name]

        for interval in tree
            region = _feature_region(genome, interval, flank; clip = false)
            isnothing(region) && continue
            feature_id, negative, region_start, region_end = region

            indices, levels, weights = _region_levels(
                find_calls_in_range(calls, max(1, region_start), region_end),
                region_start,
                region_end,
                negative,
                min_depth_bits,
                context_bits,
            )
            region_length = region_end - region_start + 1
            features[feature_id] = FeatureLevels(
                sparsevec(indices, levels, region_length),
                sparsevec(indices, weights, region_length),
            )
        end
    end

    return MethylationFrequency(min_depth_bits, context_bits, features)
end

#= Metagene profiles =#

"""
Given a vector of per-base counts, iterate the `(position, count)` pairs that
could contribute a nonzero frequency, ascending. A `SparseVector` yields its
stored entries directly (a stored zero is harmless, only wasted work); any other
vector is scanned with its zeros skipped.
"""
_nonzero_bases(counts::SparseVector) =
    zip(SparseArrays.nonzeroinds(counts), SparseArrays.nonzeros(counts))

_nonzero_bases(counts::AbstractVector) =
    ((base, count) for (base, count) in enumerate(counts) if count != 0)

"""
Given base `base` of a `2 * flank + body_length` region, find its slot in a
`2 * flank + body_bins` profile. Returns `0` for a base inside the body, which
each caller then reduces its own way; the flanks stay per base.
"""
@inline function _flank_slot(
    base::Integer,
    flank::Integer,
    body_length::Integer,
    body_bins::Integer,
)
    base <= flank && return Int(base)
    # Past the body: shift by however much the body shrank.
    base > flank + body_length && return Int(base) - Int(body_length) + Int(body_bins)
    return 0
end

"""
    gene_profile(counts, n_measurements; flank = 500, body_bins = 100)

Given one feature's per-base overlap `counts` — `flank` bp upstream, the body,
then `flank` bp downstream — reduce them to a frequency profile. Returns a
vector of length `2 * flank + body_bins`, the flanks per base and the body
interpolated onto `body_bins` evenly spaced points, every value divided by
`n_measurements`; or `nothing` when `counts` is shorter than `2 * flank + 2`.

**NOTE:** a position with no count takes `0.0 / n_measurements`, which is `NaN`
rather than `0.0` when `n_measurements` is `0`.
"""
function gene_profile(
    counts::AbstractVector,
    n_measurements::Integer;
    flank::Integer = 500,
    body_bins::Integer = 100,
)
    region_length = length(counts)
    region_length < 2 * flank + 2 && return nothing
    # Interpolating onto a single point is not defined; checked up front so a
    # feature with an empty body fails the same way as a covered one.
    body_bins >= 2 ||
        throw(ArgumentError("`body_bins` must be at least 2 (got $body_bins)"))

    body_length = region_length - 2 * flank
    unmeasured = 0.0 / n_measurements
    profile = fill(unmeasured, 2 * flank + body_bins)
    # Allocated on the first count landing in the body, and left `nothing` when
    # none does — an all-zero body interpolates to the fill above whatever its
    # length.
    body = nothing

    for (base, count) in _nonzero_bases(counts)
        frequency = Float64(count) / n_measurements
        slot = _flank_slot(base, flank, body_length, body_bins)
        if slot > 0
            profile[slot] = frequency
        else
            isnothing(body) && (body = fill(unmeasured, body_length))
            body[base-flank] = frequency
        end
    end

    if !isnothing(body)
        binned = linear_interpolation(range(0, 1; length = body_length), body).(
            range(0, 1; length = body_bins),
        )
        copyto!(profile, flank + 1, binned, 1, body_bins)
    end

    return profile
end

"""
    gene_profile(feature_levels::FeatureLevels; flank = 500, body_bins = 100,
                 weight_by_depth = true, weight_transform = identity)

Given one feature's per-base methylation, reduce it to a profile of length
`2 * flank + body_bins`. Returns a `NamedTuple` of `levels` and `weights`, or
`nothing` when the region is shorter than `2 * flank + 2`. The body is
**binned**, not interpolated: interpolating would invent levels across the long
uncovered stretches between cytosines. See `docs/weighting.md` for the weighting.

**NOTE:** `weights[i]` is zero wherever nothing was measured, and `levels[i]` is
meaningless there; a metagene needs both (see [`mean_gene_profile`](@ref)).
"""
function gene_profile(
    feature_levels::FeatureLevels;
    flank::Integer = 500,
    body_bins::Integer = 100,
    weight_by_depth::Bool = true,
    weight_transform = identity,
)
    region_length = length(feature_levels.levels)
    region_length < 2 * flank + 2 && return nothing

    n_slots = 2 * flank + body_bins
    profile_levels = zeros(Float64, n_slots)
    profile_weights = zeros(Float64, n_slots)
    weighted_sums = zeros(Float64, n_slots)

    body_length = region_length - 2 * flank
    weight_bases, base_weights = findnz(feature_levels.weights)
    level_bases, base_levels = findnz(feature_levels.levels)

    # `levels` and `weights` share a sparsity pattern in practice but need not,
    # so both stored-entry lists — each ascending — are walked in step rather
    # than indexing `levels` per base, which would binary search it every time.
    level_entry = 1
    n_level_entries = length(level_bases)

    for (entry, base) in enumerate(weight_bases)
        depth = base_weights[entry]
        depth == 0 && continue
        weight = weight_by_depth ? weight_of(depth, weight_transform) : 1.0
        weight == 0 && continue

        while level_entry <= n_level_entries && level_bases[level_entry] < base
            level_entry += 1
        end
        level =
            level_entry <= n_level_entries && level_bases[level_entry] == base ?
            Float64(base_levels[level_entry]) : 0.0

        slot = _flank_slot(base, flank, body_length, body_bins)
        slot == 0 && (slot = flank + cld((base - flank) * body_bins, body_length))
        weighted_sums[slot] += level * weight
        profile_weights[slot] += weight
    end

    for slot in eachindex(profile_levels)
        profile_weights[slot] > 0 &&
            (profile_levels[slot] = weighted_sums[slot] / profile_weights[slot])
    end

    return (levels = profile_levels, weights = profile_weights)
end

"""
    mean_gene_profile(feature_frequency; exclude = Set{String}(), flank = 500,
                      body_bins = 100)

Given per-gene coverage counts, average their [`gene_profile`](@ref)s into one
metagene profile. Returns a vector of length `2 * flank + body_bins`, all-zero
when no gene qualifies; genes listed in `exclude`, or whose stored vector is too
short for a body, are skipped.
"""
function mean_gene_profile(
    feature_frequency::FeatureFrequency;
    exclude = Set{String}(),
    flank::Integer = 500,
    body_bins::Integer = 100,
)
    accumulator = zeros(Float64, 2 * flank + body_bins)
    n_genes = 0
    for (gene_id, counts) in feature_frequency.features
        gene_id in exclude && continue
        profile = gene_profile(counts, feature_frequency.n; flank, body_bins)
        isnothing(profile) && continue
        accumulator .+= profile
        n_genes += 1
    end
    return accumulator ./ max(n_genes, 1)
end

"""
    mean_gene_profile(methylation_frequency; exclude = Set{String}(), flank = 500,
                      body_bins = 100, weight_by_depth = true,
                      weight_transform = identity)

Given per-gene methylation levels, average their [`gene_profile`](@ref)s into one
metagene profile. Returns a vector of length `2 * flank + body_bins`; genes in
`exclude`, or too short for a body, are skipped, and weighting is as in
[`gene_profile`](@ref) (see `docs/weighting.md`).

**NOTE:** each position averages only over the genes measured there, so the curve
reads as mean level among measured cytosines, not level times cytosine density;
a position no gene measured comes back as `NaN`.
"""
function mean_gene_profile(
    methylation_frequency::MethylationFrequency;
    exclude = Set{String}(),
    flank::Integer = 500,
    body_bins::Integer = 100,
    weight_by_depth::Bool = true,
    weight_transform = identity,
)
    weighted_sums = zeros(Float64, 2 * flank + body_bins)
    weight_totals = zeros(Float64, 2 * flank + body_bins)

    for (gene_id, feature_levels) in methylation_frequency.features
        gene_id in exclude && continue
        profile = gene_profile(
            feature_levels;
            flank,
            body_bins,
            weight_by_depth,
            weight_transform,
        )
        isnothing(profile) && continue

        for slot in eachindex(weighted_sums)
            # `gene_profile` already applied `weight_transform` to each base, so
            # this sum is used as-is: transforming it again would give
            # `f(Σ f(depth))`, which is not a weight on any observation.
            profile.weights[slot] == 0 && continue
            weight = weight_by_depth ? profile.weights[slot] : 1.0
            weighted_sums[slot] += profile.levels[slot] * weight
            weight_totals[slot] += weight
        end
    end

    return [
        weight_totals[slot] == 0 ? NaN : weighted_sums[slot] / weight_totals[slot] for
        slot in eachindex(weighted_sums)
    ]
end

export coverage,
    kde,
    quantiles,
    calculate_frequency,
    feature_frequency,
    FeatureFrequency,
    FeatureLevels,
    MethylationFrequency,
    DEFAULT_MIN_DEPTH,
    gene_profile,
    mean_gene_profile

end
