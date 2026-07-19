module Exploration

using Distributions
using KernelDensity
using SparseArrays
using StatsBase

using ..Reference
using ..Data

"""
    coverage(genome, feature, data::BedData; filter_zeros = false)

Return, per scaffold, how much of each feature in `genome` is covered by the
intervals in `data`.

The `feature`-type features of `genome` are intersected with `data`, and every
feature on a scaffold is scored as `covered_length / feature_length` — a value
in `[0, 1]`, where `0.0` means no overlap. Results are returned as a `Dict`
mapping scaffold name to its vector of fractions; scaffolds absent from `data`
are omitted. Set `filter_zeros = true` to drop the uncovered (`0.0`) entries.
"""
function coverage(
    genome::Genome,
    feature::Union{String,Symbol},
    data::BedData;
    filter_zeros::Bool = false,
)::Dict{String,Vector{Float64}}

    intersection = Data.intersect(genome, data, feature)
    scaffolds = Dict{String,Vector{Float64}}()
    for (scaffold_name, scaffold) in genome.scaffolds
        if !haskey(intersection, scaffold_name)
            continue
        end

        left_tree = scaffold.features
        right_tree = intersection[scaffold_name]
        iter = leftjoin(left_tree, right_tree)
        frac_coverage = Float64[]
        sizehint!(frac_coverage, length(left_tree))
        for pair in iter
            if isnothing(pair[2])
                push!(frac_coverage, 0)
            else
                length_subject = pair[1].last - pair[1].first + 1
                length_object = pair[2].last - pair[2].first + 1
                @assert length_subject >= length_object
                push!(frac_coverage, length_object / length_subject)
            end
        end

        if filter_zeros
            frac_coverage = frac_coverage |> filter(x -> x != 0)
        end

        scaffolds[scaffold_name] = frac_coverage
    end
    return scaffolds
end

"""
    kde(genome, feature, data::BedData; filter_zeros = false)

Estimate the distribution of per-feature coverage fractions for each scaffold.

Coverage is computed with [`coverage`](@ref) and a kernel density estimate is
fit to each scaffold's vector of fractions. Returns a `Dict` mapping scaffold
name to a `UnivariateKDE`, or to `nothing` when the scaffold has no fractions to
fit (for example when `filter_zeros` removed them all).
"""
function kde(
    genome::Genome,
    feature::Union{String,Symbol},
    data::BedData;
    filter_zeros::Bool = false,
)
    frac_coverage = coverage(genome, feature, data, filter_zeros = filter_zeros)
    coverage_kde = Dict{String,Union{Nothing,UnivariateKDE}}()
    for k in keys(frac_coverage)
        coverage_kde[k] =
            isempty(frac_coverage[k]) ? nothing : KernelDensity.kde(frac_coverage[k])
    end
    return coverage_kde
end

function kde(
    data::TabularData;
    filter_zeros::Bool = false,
    transform::Function = identity
)
    flat = data.table |> vec |> xs -> map(transform, xs)
    if filter_zeros
        flat = filter(x -> x != 0, flat)
    end
    flat = filter(x -> !isnan(x) && isfinite(x), flat)
    return isempty(flat) ? nothing : KernelDensity.kde(flat)
end

# Assign `value` to a 1-based bin in `1:quantiles` given the sorted quantile
# `edges` (length `quantiles + 1`). Values landing on an edge fall into the
# lower bin; anything at or above the top edge lands in the last bin.
function _quantile_bin(edges::AbstractVector, value::Real, quantiles::Int)
    for i = 1:quantiles
        value <= edges[i+1] && return i
    end
    return quantiles
end

"""
    quantiles(genome, data::TabularData; quantiles = 4, merge = mean)

Assign each sample in `data` to one of `quantiles` bins by a scalar summary of
its row.

Each matched sample's row is collapsed to a number with `merge` (default
`mean`); those values define `quantiles + 1` quantile edges, and every sample is
placed in a 1-based bin (`1` = lowest values). Returns a flat vector of
`(FeatureRecord, merged_value, quantile_index)` tuples, in the sample order of
`data`. The `FeatureRecord` is looked up in `genome` by the sample's 32-bit
metadata index; unmatched samples and unresolvable features are skipped.

# Keyword arguments
- `quantiles::Int = 4`: number of quantile bins (throws `ArgumentError` if `< 1`).
- `merge = mean`: function collapsing a row of variable values to a scalar.
"""
function quantiles(genome::Genome, data::TabularData; quantiles::Int = 4, merge = mean)
    quantiles >= 1 ||
        throw(ArgumentError("`quantiles` must be a positive integer (got $quantiles)"))

    # Collapse each matched sample's row to a scalar, remembering its metadata index.
    indices = UInt32[]
    merged = Float64[]
    for (row, sample) in enumerate(data.samples)
        sample === nothing && continue
        push!(indices, Reference.parse_index(sample[2].value))
        push!(merged, Float64(merge(data.table[row, :])))
    end

    result = Tuple{FeatureRecord,Float64,Int}[]
    isempty(merged) && return result

    # `quantiles + 1` edges spanning the observed range of merged values.
    edges = quantile(merged, range(0, 1; length = quantiles + 1))

    for (meta_idx, value) in zip(indices, merged)
        record = genome[meta_idx]
        record === nothing && continue
        push!(result, (record, value, _quantile_bin(edges, value, quantiles)))
    end

    return result
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
        error(
            "calculate_frequency supports at most $(Int(typemax(UInt16))) BedData measurements (received $n)",
        )
    end

    # Every scaffold that appears in at least one measurement. `Threads.@threads`
    # needs an indexable collection, so collect the set into a vector.
    scaffold_names = String[]
    seen = Set{String}()
    for measurement in measurements
        for name in keys(measurement.scaffolds)
            name in seen || (push!(seen, name); push!(scaffold_names, name))
        end
    end

    # Pre-populate every key so the parallel loop only overwrites existing
    # entries. Inserting new keys concurrently would race on the Dict's internal
    # structure; overwriting the value of an existing key does not.
    genome = Dict{String,SparseVector{T,Int}}(
        name => spzeros(T, 0) for name in scaffold_names
    )

    Threads.@threads for name in scaffold_names
        # Difference array: +1 where a segment starts, -1 just past its end. The
        # running total while sweeping left-to-right is the per-base frequency
        # across measurements.
        deltas = Dict{Int,Int}()
        scaffold_len = 0

        for measurement in measurements
            haskey(measurement.scaffolds, name) || continue
            tree = measurement.scaffolds[name]
            segments =
                merge ? merge_segments(tree) :
                [(Int(iv.first), Int(iv.last)) for iv in tree]
            for (s, e) in segments
                deltas[s] = get(deltas, s, 0) + 1
                deltas[e+1] = get(deltas, e + 1, 0) - 1
                scaffold_len = max(scaffold_len, e)
            end
        end

        isempty(deltas) && continue

        # Sweep the breakpoints in order, emitting a value for every covered base.
        breakpoints = sort!(collect(keys(deltas)))
        indices = Int[]
        values = T[]
        coverage = 0
        for k in eachindex(breakpoints)
            p = breakpoints[k]
            coverage += deltas[p]
            if coverage > 0 && k < length(breakpoints)
                for base = p:(breakpoints[k+1]-1)
                    push!(indices, base)
                    push!(values, T(coverage))
                end
            end
        end

        genome[name] = sparsevec(indices, values, scaffold_len)
    end

    return genome
end

"""
Per-feature, per-base coverage counts across a set of `BedData` measurements,
together with the number of measurements (`n`) they were computed from.

`features` maps a feature ID to a `SparseVector{UInt32}` of raw overlap counts —
one entry per base of the (flanked) feature, oriented in the direction of
transcription (index 1 is the feature's 5' end). Dividing a count by `n` gives
the fraction of measurements covering that base; the raw count is kept so it can
be stored exactly in 32 bits.
"""
struct FeatureFrequency
    n::Int
    features::Dict{String,SparseVector{UInt32,Int}}
end

"""
    feature_frequency(genome, feature, frequency, n; flank = 500)

Project a per-base `frequency` dictionary (as returned by
[`calculate_frequency`](@ref)) onto every `feature`-type feature of `genome`.

For each feature, the padded region `[first - flank, last + flank]` is sliced out
of its scaffold's frequency vector and re-indexed to a 1-based position within
the region. Features on the negative strand are reversed so index 1 always lands
at the feature's 5' end. The result is returned as a [`FeatureFrequency`](@ref):
a mapping from feature ID to a `SparseVector{UInt32}` of raw overlap counts, plus
the measurement count `n` so per-base frequencies can be recovered by division.
Features whose metadata ID cannot be resolved are skipped.
"""
function feature_frequency(
    genome::Genome,
    feature::Union{AbstractString,Symbol},
    frequency::AbstractDict{String,<:AbstractVector},
    n::Integer;
    flank::Integer = 500,
)
    feature_intervals = get_feature(genome, feature)
    features = Dict{String,SparseVector{UInt32,Int}}()

    for (scaffold, tree) in feature_intervals
        counts = get(frequency, scaffold, nothing)
        # Nonzero (base, count) pairs for this scaffold, ascending by position.
        nzi, nzv = counts === nothing ? (Int[], UInt32[]) : findnz(counts)

        for iv in tree
            code = iv.value
            feature_id = Reference.get_metadata_id(genome, Reference.parse_index(code))
            feature_id === nothing && continue
            negative = Reference.parse_strand(code) == get_strand('-')

            region_start = max(1, Int(iv.first) - flank)
            region_end = Int(iv.last) + flank
            region_len = region_end - region_start + 1

            # Slice of nonzero bases falling inside the padded region.
            lo = searchsortedfirst(nzi, region_start)
            hi = searchsortedlast(nzi, region_end)
            len = max(hi - lo + 1, 0)
            idxs = Vector{Int}(undef, len)
            vals = Vector{UInt32}(undef, len)
            for (j, k) in enumerate(lo:hi)
                base = nzi[k]
                # Map genomic base to a 1-based position within the region,
                # reversing for negative-strand features so index 1 stays at the
                # 5' end.
                idxs[j] = negative ? region_end - base + 1 : base - region_start + 1
                vals[j] = UInt32(nzv[k])
            end
            features[feature_id] = sparsevec(idxs, vals, region_len)
        end
    end

    return FeatureFrequency(n, features)
end

export coverage, kde, quantiles, calculate_frequency, feature_frequency, FeatureFrequency


end
