module Exploration

using Distributions
using KernelDensity
using StatsBase

using ..Reference
using ..Studies

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

    intersection = Studies.intersect(genome, data, feature)
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
    get_quantiles(genome, data::TabularData; quantiles = 4, merge = mean)

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
function get_quantiles(genome::Genome, data::TabularData; quantiles::Int = 4, merge = mean)
    quantiles >= 1 ||
        throw(ArgumentError("`quantiles` must be a positive integer (got $quantiles)"))

    # Collapse each matched sample's row to a scalar, remembering its metadata index.
    indices = UInt32[]
    merged = Float64[]
    for (row, sample) in enumerate(data.samples)
        sample === nothing && continue
        push!(indices, sample[2])
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

export coverage, kde, get_quantiles

end
