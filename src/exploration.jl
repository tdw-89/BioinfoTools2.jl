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

#= Per-feature coverage =#

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

Return, per scaffold, the mean per-base methylation level of each
`feature`-type feature of `genome`.

Each feature is scored by summing the methylation fraction of every call inside
`[first, last]` and dividing by the feature's length, giving a value in
`[0, 1]` on the same scale as the [`BedData`](@ref) method. **Bases with no call
count as zero**, so a score reflects methylation level and cytosine density
together, not level alone. Calls are summed regardless of strand and context.

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

# Fit a KDE to each scaffold's vector, mapping scaffolds with nothing to fit to
# `nothing`.
function _kde_by_scaffold(by_scaffold::Dict{String,Vector{Float64}})
    return Dict{String,Union{Nothing,UnivariateKDE}}(
        name => isempty(scores) ? nothing : KernelDensity.kde(scores) for
        (name, scores) in by_scaffold
    )
end

"""
    kde(data::BedData, feature; filter_zeros = false)

Estimate the distribution of per-feature coverage fractions for each scaffold.

Coverage is computed with [`coverage`](@ref) and a kernel density estimate is
fit to each scaffold's vector of fractions. Returns a `Dict` mapping scaffold
name to a `UnivariateKDE`, or to `nothing` when the scaffold has no fractions to
fit (for example when `filter_zeros` removed them all).
"""
kde(data::BedData, feature::Union{AbstractString,Symbol}; filter_zeros::Bool = false) =
    _kde_by_scaffold(coverage(data, feature; filter_zeros = filter_zeros))

"""
    kde(genome, data::MethylationData, feature; filter_zeros = false)

Estimate the distribution of per-feature mean methylation levels for each
scaffold.

Levels are computed with [`coverage`](@ref) — fractions in `[0, 1]`, with
uncovered bases counted as zero — and a kernel density estimate is fit to each
scaffold's vector. Returns a `Dict` mapping scaffold name to a `UnivariateKDE`,
or to `nothing` when the scaffold has nothing to fit.
"""
kde(
    genome::Genome,
    data::MethylationData,
    feature::Union{AbstractString,Symbol};
    filter_zeros::Bool = false,
) = _kde_by_scaffold(coverage(genome, data, feature; filter_zeros = filter_zeros))

"""
    kde(data::TabularData; filter_zeros = false, transform = identity)

Estimate the distribution of every value in `data`'s table, after applying
`transform`. Non-finite values are dropped, as are zeros when
`filter_zeros = true`. Returns `nothing` when nothing is left to fit.
"""
function kde(data::TabularData; filter_zeros::Bool = false, transform::Function = identity)
    flat = [transform(value) for value in data.table]
    filter!(value -> isfinite(value) && !(filter_zeros && iszero(value)), flat)
    return isempty(flat) ? nothing : KernelDensity.kde(flat)
end

#= Quantile binning =#

# Assign `value` to a 1-based bin in `1:n_bins` given the sorted quantile `edges`
# (length `n_bins + 1`). Values landing on an edge fall into the lower bin;
# anything at or above the top edge lands in the last bin.
_quantile_bin(edges::AbstractVector, value::Real, n_bins::Int) =
    min(searchsortedfirst(view(edges, 2:lastindex(edges)), value), n_bins)

# Validate the arguments the rank-based `quantiles` methods share. `noun` names
# what `ranking` selects, for the error message.
function _check_ranking(n_bins::Int, ranking::Vector{String}, noun::String)
    n_bins >= 1 ||
        throw(ArgumentError("`quantiles` must be a positive integer (got $n_bins)"))
    isempty(ranking) && throw(ArgumentError("`ranking` must name at least one $noun"))
    return nothing
end

# Sort `ranking_keys` — tuples of ranking values with the original row index
# last — and return `(row, bin)` pairs in ranked order, cutting the ranking into
# `n_bins` equal-sized (up to an off-by-one) bins.
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

Assign each sample in `data` to one of `quantiles` bins by a scalar summary of
its row.

Each matched sample's row is collapsed to a number with `merge` (default
`mean`); those values define `quantiles + 1` quantile edges, and every sample is
placed in a 1-based bin (`1` = lowest values). Returns a flat vector of
`(FeatureRecord, merged_value, quantile_index)` tuples, in the sample order of
`data`. The `FeatureRecord` is looked up in `data.genome` by the sample's 32-bit
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

Rank samples in `data` by the named `ranking` variables — the first entry
primary, each subsequent entry breaking ties in the one before it, and any
remaining tie broken by original sample order — then cut that ranking into
`quantiles` equal-sized (up to an off-by-one) bins.

Returns a flat vector of `(FeatureRecord, quantile_index)` tuples, in ranked
order. The `FeatureRecord` is looked up in `data.genome` by the sample's 32-bit
metadata index; unmatched samples and unresolvable features are skipped, but
(like every matched sample) still occupy a rank position and so still affect
bin sizing.
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

Rank the rows of a paralog-pair table (shaped like `GeneFamily`'s constructor
input, or `rbh`'s output) by the named `ranking` columns — the first entry
primary, each subsequent entry breaking ties in the one before it, and any
remaining tie broken by original row order — then cut that ranking into
`quantiles` equal-sized (up to an off-by-one) bins.

Returns `pairs` with a `"quantile"` column appended (1-based bin number); row
order is unchanged.
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
    calculate_frequency(measurements; merge = true)

Count how many of the `BedData` `measurements` cover each base of the genome,
returning a `Dict` mapping scaffold name to a `SparseVector` of per-base
frequencies (length = the largest interval end seen on that scaffold).

When `merge` is `true` (the default), overlapping intervals *within a single
measurement* are merged first (via [`merge_segments`](@ref)), so each
measurement contributes at most 1 to a given base and the maximum possible value
is the number of measurements. Set `merge = false` to skip this step when the
intervals are already disjoint (e.g. ChIP-seq peak calls), in which case any
within-measurement overlaps will stack.

The element type is chosen to fit the measurement count: `UInt8` for up to 255
measurements and `UInt16` for up to 65535. More measurements raise an error.
"""
function calculate_frequency(measurements::Vector{BedData}; merge::Bool = true)
    n_measurements = length(measurements)
    T = if n_measurements <= typemax(UInt8)
        UInt8
    elseif n_measurements <= typemax(UInt16)
        UInt16
    else
        error(
            "calculate_frequency supports at most $(Int(typemax(UInt16))) BedData measurements (received $n_measurements)",
        )
    end

    # Every scaffold that appears in at least one measurement. `Threads.@threads`
    # needs an indexable collection, so collect the set into a vector.
    scaffold_names = String[]
    seen = Set{String}()
    for measurement in measurements, name in keys(measurement.scaffolds)
        name in seen || (push!(seen, name); push!(scaffold_names, name))
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
        events = Tuple{Int,Int}[]
        scaffold_length = 0

        for measurement in measurements
            tree = get(measurement.scaffolds, name, nothing)
            isnothing(tree) && continue
            segments =
                merge ? merge_segments(tree) :
                [(Int(iv.first), Int(iv.last)) for iv in tree]
            for (start_pos, end_pos) in segments
                push!(events, (start_pos, 1))
                push!(events, (end_pos + 1, -1))
                scaffold_length = max(scaffold_length, end_pos)
            end
        end

        isempty(events) && continue
        sort!(events)

        # Sweep the breakpoints in order, emitting a value for every covered base.
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

        genome[name] = sparsevec(indices, counts, scaffold_length)
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

Project a per-base `frequency` dictionary (as returned by
[`calculate_frequency`](@ref)) onto every `feature`-type feature of `genome`.

For each feature, the padded region `[first - flank, last + flank]` is sliced out
of its scaffold's frequency vector and re-indexed to a 1-based position within
the region. Features on the negative strand are reversed so index 1 always lands
at the feature's 5' end. The result is returned as a [`FeatureFrequency`](@ref),
carrying the measurement count `n` so per-base frequencies can be recovered by
division. Features whose metadata ID cannot be resolved are skipped.
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
            code = interval.value
            feature_id = Reference.get_metadata_id(genome, Reference.parse_index(code))
            isnothing(feature_id) && continue
            negative = Reference.parse_strand(code) == get_strand('-')

            region_start = max(1, Int(interval.first) - flank)
            region_end = Int(interval.last) + flank

            # Slice of nonzero bases falling inside the padded region.
            first_entry = searchsortedfirst(base_indices, region_start)
            last_entry = searchsortedlast(base_indices, region_end)
            n_entries = max(last_entry - first_entry + 1, 0)
            indices = Vector{Int}(undef, n_entries)
            counts_in_region = Vector{UInt32}(undef, n_entries)
            for (slot, entry) in enumerate(first_entry:last_entry)
                base = base_indices[entry]
                # Map genomic base to a 1-based position within the region,
                # reversing for negative-strand features so index 1 stays at the
                # 5' end.
                indices[slot] = negative ? region_end - base + 1 : base - region_start + 1
                counts_in_region[slot] = UInt32(base_counts[entry])
            end
            features[feature_id] =
                sparsevec(indices, counts_in_region, region_end - region_start + 1)
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

Turn a read depth into the statistical weight the profile functions give it, by
applying `weight_transform` (`identity` for linear weighting, `log`, `sqrt`, …).
The depth is converted to `Float64` first, since depths are stored unsigned.

Throws an `ArgumentError` for a weight that is negative or not finite. Note that
a transform may legitimately return `0` — `log(1) == 0` — which drops that
observation; prefer `log1p` over `log` unless that is what you want.
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
    context_bits = isnothing(context) ? nothing : UInt8(context)
    features = Dict{String,FeatureLevels}()

    for (scaffold_name, tree) in get_feature(genome, feature)
        haskey(data, scaffold_name) || continue
        calls = data[scaffold_name]

        for interval in tree
            code = interval.value
            feature_id = Reference.get_metadata_id(genome, Reference.parse_index(code))
            isnothing(feature_id) && continue
            negative = Reference.parse_strand(code) == get_strand('-')

            region_start = Int(interval.first) - flank
            region_end = Int(interval.last) + flank

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
                    isnothing(context_bits) || get_context(call) == context_bits || continue
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

            region_length = region_end - region_start + 1
            features[feature_id] = FeatureLevels(
                sparsevec(indices, level_values, region_length),
                sparsevec(indices, weight_values, region_length),
            )
        end
    end

    return MethylationFrequency(min_depth_bits, context_bits, features)
end

#= Metagene profiles =#

"""
Iterate the `(position, count)` pairs of `counts` that could contribute a nonzero
frequency, in ascending 1-based position. A `SparseVector` yields its stored
entries directly (a stored zero is harmless, only wasted work); any other vector
is scanned with its zeros skipped.
"""
_nonzero_bases(counts::SparseVector) =
    zip(SparseArrays.nonzeroinds(counts), SparseArrays.nonzeros(counts))

_nonzero_bases(counts::AbstractVector) =
    ((base, count) for (base, count) in enumerate(counts) if count != 0)

# Position in a `2 * flank + body_bins` profile that base `base` of a
# `2 * flank + body_length` region falls in. Flanks stay per base; a body base
# is mapped by `body_slot`, which either interpolates or bins it.
@inline function _profile_slot(
    base::Int,
    flank::Int,
    body_length::Int,
    body_slot::F,
) where {F}
    base <= flank && return base
    # Past the body: shift by however much the body shrank.
    base > flank + body_length && return base - body_length
    return -1
end

"""
    gene_profile(counts, n_measurements; flank = 500, body_bins = 100)

Reduce one feature's per-base overlap `counts` — `flank` bp upstream, the feature
body, then `flank` bp downstream — to a frequency profile of length
`2 * flank + body_bins`. The flanks are kept per base while the body is
interpolated onto `body_bins` evenly spaced points, and every value is divided by
`n_measurements` to give a frequency. Returns `nothing` when `counts` is shorter
than `2 * flank + 2` (no room for a body of at least two bases).

`counts` is typically a `SparseVector` from [`FeatureFrequency`](@ref) holding a
few thousand nonzeros in a region tens of kilobases long, and most features carry
none at all, so only the stored entries are visited and the body is materialised
densely only when one falls inside it. Positions with no count take
`0.0 / n_measurements`, which is `NaN` rather than `0.0` when `n_measurements`
is `0`.
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
        if base <= flank
            profile[base] = frequency
        elseif base > flank + body_length
            # Past the body: shift by however much the body shrank.
            profile[base-body_length+body_bins] = frequency
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

        # Flanks stay per base; body bases fold into one of `body_bins` bins.
        slot = if base <= flank
            base
        elseif base > flank + body_length
            # Past the body: shift by however much the body shrank.
            base - body_length + body_bins
        else
            flank + cld((base - flank) * body_bins, body_length)
        end

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
same scheme. See [`weight_of`](@ref) and `docs/weighting.md`.

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
