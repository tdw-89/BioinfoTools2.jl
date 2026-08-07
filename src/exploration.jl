module Exploration

# Scoped, not a blanket `using DataFrames` — that would also pull in
# `DataFrames.leftjoin`, colliding with `Data.leftjoin` used by `coverage`.
using DataFrames: DataFrame, nrow
using Distributions
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

"""
    coverage(data::BedData, feature; filter_zeros = false)

Return, per scaffold, how much of each feature in `data.genome` is covered by the
intervals in `data`.

The `feature`-type features of `genome` are intersected with `data`, and every
feature on a scaffold is scored as `covered_length / feature_length` — a value
in `[0, 1]`, where `0.0` means no overlap. Results are returned as a `Dict`
mapping scaffold name to its vector of fractions; scaffolds absent from `data`
are omitted. Set `filter_zeros = true` to drop the uncovered (`0.0`) entries.
"""
function coverage(
    data::BedData,
    feature::Union{String,Symbol};
    filter_zeros::Bool = false,
)::Dict{String,Vector{Float64}}

    intersection = Data.intersect(data.genome, data, feature)
    scaffolds = Dict{String,Vector{Float64}}()
    for (scaffold_name, scaffold) in data.genome.scaffolds
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
    coverage(genome, data::MethylationData, feature; filter_zeros = false)

Return, per scaffold, the mean per-base methylation level of each
`feature`-type feature of `genome`.

Each feature is scored by summing the methylation fraction of every call inside
`[first, last]` and dividing by the feature's length, giving a value in
`[0, 1]` on the same scale as the [`BedData`](@ref) method. **Bases with no call
count as zero**: cytosines are sparse, so a score reflects methylation level and
cytosine density together, not level alone. Calls are summed regardless of
strand and context.

Results are returned as a `Dict` mapping scaffold name to its vector of scores;
scaffolds absent from `data` are omitted. Set `filter_zeros = true` to drop
features with no covered cytosine at all.
"""
function coverage(
    genome::Genome,
    data::MethylationData,
    feature::Union{AbstractString,Symbol};
    filter_zeros::Bool = false,
)::Dict{String,Vector{Float64}}

    feature_intervals = get_feature(genome, feature)
    scaffolds = Dict{String,Vector{Float64}}()
    for (scaffold_name, tree) in feature_intervals
        haskey(data, scaffold_name) || continue
        calls = data[scaffold_name]

        mean_levels = Float64[]
        sizehint!(mean_levels, length(tree))
        for interval in tree
            feature_length = Int(interval.last) - Int(interval.first) + 1
            total_fraction = 0.0
            for call in find_calls_in_range(calls, interval.first, interval.last)
                fraction = meth_fraction(call)
                # `meth_fraction` is NaN for a site with no coverage at all.
                isnan(fraction) || (total_fraction += fraction)
            end
            push!(mean_levels, total_fraction / feature_length)
        end

        if filter_zeros
            mean_levels = filter(!iszero, mean_levels)
        end
        scaffolds[scaffold_name] = mean_levels
    end
    return scaffolds
end

"""
    kde(data::BedData, feature; filter_zeros = false)

Estimate the distribution of per-feature coverage fractions for each scaffold.

Coverage is computed with [`coverage`](@ref) and a kernel density estimate is
fit to each scaffold's vector of fractions. Returns a `Dict` mapping scaffold
name to a `UnivariateKDE`, or to `nothing` when the scaffold has no fractions to
fit (for example when `filter_zeros` removed them all).
"""
function kde(data::BedData, feature::Union{String,Symbol}; filter_zeros::Bool = false)
    frac_coverage = coverage(data, feature, filter_zeros = filter_zeros)
    coverage_kde = Dict{String,Union{Nothing,UnivariateKDE}}()
    for k in keys(frac_coverage)
        coverage_kde[k] =
            isempty(frac_coverage[k]) ? nothing : KernelDensity.kde(frac_coverage[k])
    end
    return coverage_kde
end

"""
    kde(genome, data::MethylationData, feature; filter_zeros = false)

Estimate the distribution of per-feature mean methylation levels for each
scaffold.

Levels are computed with [`coverage`](@ref) — fractions in `[0, 1]`, with
uncovered bases counted as zero — and a kernel density estimate is fit to each
scaffold's vector. Returns a `Dict` mapping scaffold name to a `UnivariateKDE`,
or to `nothing` when the scaffold has nothing to fit.
"""
function kde(
    genome::Genome,
    data::MethylationData,
    feature::Union{AbstractString,Symbol};
    filter_zeros::Bool = false,
)
    mean_levels = coverage(genome, data, feature; filter_zeros = filter_zeros)
    level_kde = Dict{String,Union{Nothing,UnivariateKDE}}()
    for (scaffold_name, levels) in mean_levels
        level_kde[scaffold_name] = isempty(levels) ? nothing : KernelDensity.kde(levels)
    end
    return level_kde
end

function kde(data::TabularData; filter_zeros::Bool = false, transform::Function = identity)
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
function quantiles(data::TabularData; quantiles::Int = 4, merge = mean)
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
        record = data.genome[meta_idx]
        record === nothing && continue
        push!(result, (record, value, _quantile_bin(edges, value, quantiles)))
    end

    return result
end

"""
    quantiles(data::TabularData, ranking::Vector{String}; quantiles::Int = 4)

Rank samples in `data` by the named `ranking` variables — the first entry
primary, each subsequent entry breaking ties in the one before it, and any
remaining tie broken by original sample order — then cut that ranking into
`quantiles` equal-sized (up to an off-by-one) bins.

Returns a flat vector of `(FeatureRecord, quantile_index)` tuples, in ranked
order. The `FeatureRecord` is looked up in `genome` by the sample's 32-bit
metadata index; unmatched samples and unresolvable features are skipped, but
(like every matched sample) still occupy a rank position and so still affect
bin sizing.
"""
function quantiles(data::TabularData, ranking::Vector{String}; quantiles::Int = 4)
    quantiles >= 1 ||
        throw(ArgumentError("`quantiles` must be a positive integer (got $quantiles)"))
    isempty(ranking) && throw(ArgumentError("`ranking` must name at least one variable"))

    col_indices = map(ranking) do name
        idx = findfirst(==(name), data.variables)
        idx === nothing && throw(ArgumentError("Unknown ranking variable: \"$name\""))
        idx
    end

    # Matched samples only; the row index itself is the final tie-break.
    rows = findall(!isnothing, data.samples)
    n = length(rows)
    n == 0 && return Tuple{FeatureRecord,Int}[]

    keyed = sort([(Tuple(data.table[row, c] for c in col_indices)..., row) for row in rows])

    result = Tuple{FeatureRecord,Int}[]
    sizehint!(result, n)
    for (rank_pos, key) in enumerate(keyed)
        row = key[end]
        meta_idx = Reference.parse_index(data.samples[row][2].value)
        record = data.genome[meta_idx]
        record === nothing && continue
        push!(result, (record, cld(rank_pos * quantiles, n)))
    end
    return result
end

"""
    quantiles(pairs::DataFrame, ranking::Vector{String}; quantiles::Int = 4) -> DataFrame

Rank the rows of a paralog-pair table (shaped like `GeneFamily`'s constructor
input, or `rbh`'s output) by the named `ranking` columns — the first entry
primary, each subsequent entry breaking ties in the one before it, and any
remaining tie broken by original row order — then cut that ranking into
`quantiles` equal-sized (up to an off-by-one) bins.

Returns `pairs` with a `"quantile"` column appended (1-based bin number); row
order is unchanged.
"""
function quantiles(pairs::DataFrame, ranking::Vector{String}; quantiles::Int = 4)
    quantiles >= 1 ||
        throw(ArgumentError("`quantiles` must be a positive integer (got $quantiles)"))
    isempty(ranking) && throw(ArgumentError("`ranking` must name at least one column"))
    for name in ranking
        name in names(pairs) || throw(ArgumentError("Unknown ranking column: \"$name\""))
    end

    n = nrow(pairs)
    bins = Vector{Int}(undef, n)
    if n > 0
        order = sort([(Tuple(pairs[row, name] for name in ranking)..., row) for row = 1:n])
        for (rank_pos, key) in enumerate(order)
            bins[key[end]] = cld(rank_pos * quantiles, n)
        end
    end

    result = copy(pairs)
    result.quantile = bins
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
    genome =
        Dict{String,SparseVector{T,Int}}(name => spzeros(T, 0) for name in scaffold_names)

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

"""
Minimum read depth a call must have to enter a methylation profile. Shallow
sites carry a level that is mostly sampling noise — a 1-read site reads as 0% or
100% — so they are dropped rather than weighted down.
"""
const DEFAULT_MIN_DEPTH = UInt32(5)

"""
    weight_of(depth, weight_transform)

Turn a read depth into the statistical weight the profile functions give it, by
applying `weight_transform` (`identity` for linear weighting, `log`, `sqrt`, …).

The depth is converted to `Float64` before the transform sees it: depths are
stored unsigned, where an expression as ordinary as `-depth` would wrap to a
huge positive number instead of going negative.

Throws an `ArgumentError` for a weight that is negative or not finite, which
would otherwise silently corrupt a weighted mean. Note that a transform may
legitimately return `0` — `log(1) == 0` — which drops that observation; prefer
`log1p` over `log` unless that is what you want.
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

The methylation counterpart of [`FeatureFrequency`](@ref). It is a separate type
because the two carry different quantities: `FeatureFrequency` holds overlap
*counts* to be divided by a measurement count, whereas this holds a *level* per
base plus the depth it rests on, and needs to distinguish an unmeasured base
from one measured at zero.
"""
struct MethylationFrequency
    min_depth::UInt32
    context::Union{Nothing,UInt8}
    features::Dict{String,FeatureLevels}
end

"""
    feature_frequency(genome, feature, data::MethylationData; flank = 500,
                      min_depth = DEFAULT_MIN_DEPTH, context = CTX_CPG)

Project the per-base methylation levels in `data` onto every `feature`-type
feature of `genome`, returning a [`MethylationFrequency`](@ref).

No [`calculate_frequency`](@ref) step is needed: a `MethylationData` already
holds one aggregated call per base. Each feature's padded region
`[first - flank, last + flank]` is sliced out and re-indexed to a 1-based
position within the region, reversed for negative-strand features so index 1 is
always the 5' end. Unlike the [`BedData`](@ref) pipeline the region is *not*
clipped at the scaffold start, so index `flank + 1` is the feature's first base
for every feature; a clipped region simply holds no data there.

Calls are dropped unless their depth is at least `min_depth` and their context
matches `context` (pass `context = nothing` to keep every context). Where two
calls share a base — the same position on opposite strands, say — their levels
are combined in proportion to their depths and their depths are summed.
Features whose metadata ID cannot be resolved are skipped.
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
    context_bits = context === nothing ? nothing : UInt8(context)

    feature_intervals = get_feature(genome, feature)
    features = Dict{String,FeatureLevels}()

    for (scaffold_name, tree) in feature_intervals
        haskey(data, scaffold_name) || continue
        calls = data[scaffold_name]

        for interval in tree
            code = interval.value
            feature_id = Reference.get_metadata_id(genome, Reference.parse_index(code))
            feature_id === nothing && continue
            negative = Reference.parse_strand(code) == get_strand('-')

            region_start = Int(interval.first) - flank
            region_end = Int(interval.last) + flank
            region_length = region_end - region_start + 1

            region_calls = find_calls_in_range(calls, max(1, region_start), region_end)
            indices = Int[]
            level_values = Float32[]
            weight_values = UInt32[]

            # Calls are position-sorted, so the (rare) several calls sharing a
            # base form one contiguous run.
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
                    depth < min_depth_bits && continue
                    context_bits === nothing ||
                        get_context(call) == context_bits ||
                        continue
                    fraction = meth_fraction(call)
                    isnan(fraction) && continue
                    # Always linear in depth: this pools reads at one cytosine
                    # back into that cytosine's own methylation fraction, which
                    # only comes out right weighted by read count.
                    weighted_fraction += fraction * Int(depth)
                    depth_total += Int(depth)
                end
                depth_total == 0 && continue

                push!(
                    indices,
                    negative ? region_end - Int(position) + 1 :
                    Int(position) - region_start + 1,
                )
                push!(level_values, Float32(weighted_fraction / depth_total))
                push!(weight_values, UInt32(min(depth_total, Int(typemax(UInt32)))))
            end

            features[feature_id] = FeatureLevels(
                sparsevec(indices, level_values, region_length),
                sparsevec(indices, weight_values, region_length),
            )
        end
    end

    return MethylationFrequency(min_depth_bits, context_bits, features)
end

"""
    gene_profile(counts, n_measurements; flank = 500, body_bins = 100)

Reduce one feature's per-base overlap `counts` — `flank` bp upstream, the feature
body, then `flank` bp downstream — to a frequency profile of length
`2 * flank + body_bins`. The flanks are kept per base while the body is
interpolated onto `body_bins` evenly spaced points, and every value is divided by
`n_measurements` to give a frequency. Returns `nothing` when `counts` is shorter
than `2 * flank + 2` (no room for a body of at least two bases).
"""
function gene_profile(
    counts::AbstractVector,
    n_measurements::Integer;
    flank::Integer = 500,
    body_bins::Integer = 100,
)
    length(counts) < 2 * flank + 2 && return nothing
    frequency = Vector{Float64}(counts) ./ n_measurements
    body = frequency[(flank+1):(end-flank)]
    body_binned = linear_interpolation(range(0, 1; length = length(body)), body).(
        range(0, 1; length = body_bins),
    )
    return vcat(frequency[1:flank], body_binned, frequency[(end-flank+1):end])
end

"""
    gene_profile(feature_levels::FeatureLevels; flank = 500, body_bins = 100,
                 weight_by_depth = true, weight_transform = identity)

Reduce one feature's per-base methylation — `flank` bp upstream, the feature
body, then `flank` bp downstream — to a profile of length
`2 * flank + body_bins`. Returns a `NamedTuple` of two vectors, `levels` and
`weights`, or `nothing` when the region is shorter than `2 * flank + 2`.

`weights[i]` is the total weight behind position `i` and is **zero wherever
nothing was measured**; `levels[i]` is the weighted mean fraction there, and is
meaningless when the weight is zero. Both are needed because a metagene averages
only over the features that carry a call at a position — see
[`mean_gene_profile`](@ref).

The flanks are kept per base. The body is **binned**, not interpolated: each
body base falls in one of `body_bins` equal bins and a bin's level is the
weighted mean of the calls landing in it. Interpolating instead, as the
[`BedData`](@ref) method does, would invent levels for the long uncovered
stretches between cytosines.

A base's weight is `weight_transform(depth)` — see [`weight_of`](@ref) — so
`identity` (the default) weights linearly by read depth and `log`/`sqrt`
compress the advantage of deeply covered sites. With `weight_by_depth = false`
every measured base counts equally and `weight_transform` is ignored.
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

    profile_levels = zeros(Float64, 2 * flank + body_bins)
    profile_weights = zeros(Float64, 2 * flank + body_bins)
    weighted_sums = zeros(Float64, 2 * flank + body_bins)

    body_length = region_length - 2 * flank
    nonzero_indices, nonzero_weights = findnz(feature_levels.weights)

    for (entry, base) in enumerate(nonzero_indices)
        depth = nonzero_weights[entry]
        depth == 0 && continue
        weight = weight_by_depth ? weight_of(depth, weight_transform) : 1.0
        weight == 0 && continue

        # Flanks stay per base; body bases fold into one of `body_bins` bins.
        slot = if base <= flank
            base
        elseif base > flank + body_length
            # Past the body: shift by however much the body shrank.
            base - body_length + body_bins
        else
            flank + cld((base - flank) * body_bins, body_length)
        end

        weighted_sums[slot] += Float64(feature_levels.levels[base]) * weight
        profile_weights[slot] += weight
    end

    for slot in eachindex(profile_levels)
        profile_weights[slot] > 0 &&
            (profile_levels[slot] = weighted_sums[slot] / profile_weights[slot])
    end

    return (levels = profile_levels, weights = profile_weights)
end

"""
    mean_gene_profile(feature_frequency; exclude = Set{String}(), flank = 500, body_bins = 100)

Average the per-gene [`gene_profile`](@ref)s in `feature_frequency` into a single
metagene profile of length `2 * flank + body_bins`. Genes listed in `exclude`, or
whose stored vector is too short for a body, are skipped. Returns an all-zero
profile when no gene qualifies.
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
        profile === nothing && continue
        accumulator .+= profile
        n_genes += 1
    end
    return accumulator ./ max(n_genes, 1)
end

"""
    mean_gene_profile(methylation_frequency; exclude = Set{String}(), flank = 500,
                      body_bins = 100, weight_by_depth = true,
                      weight_transform = identity)

Average the per-gene [`gene_profile`](@ref)s in `methylation_frequency` into a
single metagene profile of length `2 * flank + body_bins`.

**Each position averages only over the genes measured there** — a gene with no
call at a flank base, or no call anywhere in a body bin, does not contribute to
that position at all rather than contributing a zero. The profile therefore
reads as the mean methylation level among measured cytosines, not as level times
cytosine density. Positions no gene measured come back as `NaN`.

With `weight_by_depth = true` (the default) a position is weighted by the read
depth behind it, both when a gene's own calls are pooled into a bin and when
genes are averaged together; `false` weights every measured base and every
contributing gene equally. `weight_transform` reshapes depth into weight —
`identity` for linear, `log`/`sqrt` to compress the advantage of deeply covered
sites — and is applied once per base, so both levels of averaging inherit the
same scheme. See [`weight_of`](@ref).

Genes listed in `exclude`, or whose region is too short for a body, are skipped.
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
        profile === nothing && continue

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
