"""
Per-base empirical mappability from UniMap.jl, as runs of identical counts.

A positional signal like [`Methylation`](@ref): nothing here touches a `Genome`.
A UniMap `mappability.zarr` is converted once into a cache of Arrow files — one
per scaffold per channel `(L, F, mapq_cutoff)` — by [`convert_unimap`](@ref)
(load `Zarr` to enable it), and read back lazily by [`read_mappability`](@ref).
See `docs/mappability.md` for the layout and the reasons for it.
"""
module Mappability

using Arrow
using StructArrays
using TOML

using ..Methylation: MethylationData, AggregatedCall, sanitize_name, _unique_stem!

#= Types =#

"""
A run of consecutive bases with identical UniMap counts, exactly 16 bytes:

| field | type | meaning |
|---|---|---|
| `start`, `stop` | `UInt32` | 1-based closed span on the scaffold |
| `n_simulated` | `UInt16` | simulated fragments covering each base |
| `n_correct` | `UInt16` | of those, placed back correctly at MAPQ ≥ cutoff |
| `n_ectopic` | `UInt32` | misplaced fragments landing on each base at MAPQ ≥ cutoff |

A base absent from every run had all three counts zero (unsimulated, or N).
"""
struct MappabilityRun
    start::UInt32
    stop::UInt32
    n_simulated::UInt16
    n_correct::UInt16
    n_ectopic::UInt32
end

"""The grid point `(L, F, mapq_cutoff)` a set of runs belongs to."""
const MappabilityChannel = @NamedTuple{L::Int, F::Int, mapq_cutoff::Int}

"""
One channel's runs, per scaffold, each a `StructArray{MappabilityRun}` sorted
by `start` and non-overlapping (so `stop` is sorted too), plus the `channel` and
the `pseudocount` [`log2_ratio`](@ref) uses.
"""
struct MappabilityData
    scaffolds::Dict{String,StructArray{MappabilityRun}}
    channel::MappabilityChannel
    pseudocount::Float64
end

Base.getindex(data::MappabilityData, scaffold::AbstractString) =
    data.scaffolds[String(scaffold)]
Base.haskey(data::MappabilityData, scaffold::AbstractString) =
    haskey(data.scaffolds, String(scaffold))
Base.keys(data::MappabilityData) = keys(data.scaffolds)
Base.length(data::MappabilityData) = length(data.scaffolds)

"""Total number of runs across every scaffold."""
n_runs(data::MappabilityData) = sum(length, values(data.scaffolds); init = 0)

#= Metrics =#

"""
Given a run, the fraction of its simulated fragments placed back correctly.
Returns `NaN32` when nothing was simulated there.
"""
map_back(run) = run.n_simulated == 0 ? NaN32 : Float32(run.n_correct / run.n_simulated)

"""
Given a run, misplaced fragments landing on it per fragment simulated there.
Returns `NaN32` when nothing was simulated there.
"""
ectopic_rate(run) = run.n_simulated == 0 ? NaN32 : Float32(run.n_ectopic / run.n_simulated)

"""
Given a run, `log2((n_correct + pseudocount) / (n_ectopic + pseudocount))`, the
log2 odds that a kept fragment there is native (UniMap's definition). Returns
`NaN32` when both counts are zero.
"""
log2_ratio(run, pseudocount::Real) =
    run.n_correct == 0 && run.n_ectopic == 0 ? NaN32 :
    Float32(log2((run.n_correct + pseudocount) / (run.n_ectopic + pseudocount)))

"""
Given a metric name — `:map_back`, `:ectopic_rate` or `:log2_ratio` — and the
data's `pseudocount`, return a function from a run to that metric. Build it
once, outside a loop: the `Symbol` dispatch is not free.
"""
function metric_function(metric::Symbol, pseudocount::Real)
    metric === :map_back && return map_back
    metric === :ectopic_rate && return ectopic_rate
    metric === :log2_ratio && return run -> log2_ratio(run, pseudocount)
    throw(
        ArgumentError(
            "unknown metric :$metric (use :map_back, :ectopic_rate or :log2_ratio)",
        ),
    )
end

#= Queries =#

"""
    runs_in_range(runs, start_pos, stop_pos)

Given one scaffold's runs, find those overlapping `[start_pos, stop_pos]`.
Returns a view (nothing copied), empty when none overlap or `stop_pos < start_pos`.
Binary searches `stop` for the first and `start` for the last, `O(log n)`.
"""
function runs_in_range(
    runs::StructArray{MappabilityRun},
    start_pos::Integer,
    stop_pos::Integer,
)
    stop_pos < start_pos && return view(runs, 1:0)
    return view(
        runs,
        searchsortedfirst(runs.stop, start_pos):searchsortedlast(runs.start, stop_pos),
    )
end

"""
    runs_in_range(data, scaffold, start_pos, stop_pos)

Range query against one scaffold of a [`MappabilityData`](@ref). Returns
`nothing` when `scaffold` isn't present.
"""
function runs_in_range(
    data::MappabilityData,
    scaffold::AbstractString,
    start_pos::Integer,
    stop_pos::Integer,
)
    runs = get(data.scaffolds, String(scaffold), nothing)
    runs === nothing && return nothing
    return runs_in_range(runs, start_pos, stop_pos)
end

"""
    mappability_at(data, scaffold, positions; metric = :map_back)

Given 1-based `positions` on `scaffold`, look up `metric` at each. Returns a
`Vector{Float32}`, `NaN` at a position no run covers (or on a missing scaffold).
"""
function mappability_at(
    data::MappabilityData,
    scaffold::AbstractString,
    positions::AbstractVector{<:Integer};
    metric::Symbol = :map_back,
)
    value_of = metric_function(metric, data.pseudocount)
    values = fill(NaN32, length(positions))
    haskey(data, scaffold) || return values
    runs = data[scaffold]
    for (slot, position) in enumerate(positions)
        entry = searchsortedlast(runs.start, position)
        entry > 0 && runs.stop[entry] >= position && (values[slot] = value_of(runs[entry]))
    end
    return values
end

"""
    per_base(data, scaffold, start_pos, stop_pos; metric = :map_back)

Given a span of `scaffold`, expand its runs to one `metric` value per base.
Returns a `Vector{Float32}` of length `stop_pos - start_pos + 1`, `NaN` where
no run covers a base.
"""
function per_base(
    data::MappabilityData,
    scaffold::AbstractString,
    start_pos::Integer,
    stop_pos::Integer;
    metric::Symbol = :map_back,
)
    value_of = metric_function(metric, data.pseudocount)
    values = fill(NaN32, max(0, stop_pos - start_pos + 1))
    runs = runs_in_range(data, scaffold, start_pos, stop_pos)
    isnothing(runs) && return values
    for run in runs
        first_base = max(Int(run.start), start_pos)
        last_base = min(Int(run.stop), stop_pos)
        values[(first_base-start_pos+1):(last_base-start_pos+1)] .= value_of(run)
    end
    return values
end

#= Building runs =#

"""
    append_runs!(runs, first_base, n_simulated, n_correct, n_ectopic)

Given per-base counts for consecutive bases starting at 1-based `first_base`,
extend `runs` with them: a base equal to, and adjacent to, the last run lengthens
it; an all-zero base is skipped. Returns `runs`. Feed bases in ascending order.
"""
function append_runs!(
    runs::Vector{MappabilityRun},
    first_base::Integer,
    n_simulated::AbstractVector,
    n_correct::AbstractVector,
    n_ectopic::AbstractVector,
)
    for offset in eachindex(n_simulated, n_correct, n_ectopic)
        simulated, correct, ectopic =
            n_simulated[offset], n_correct[offset], n_ectopic[offset]
        simulated == 0 && ectopic == 0 && continue
        base = UInt32(first_base + offset - firstindex(n_simulated))
        if !isempty(runs)
            last_run = runs[end]
            if last_run.stop + 1 == base &&
               last_run.n_simulated == simulated &&
               last_run.n_correct == correct &&
               last_run.n_ectopic == ectopic
                runs[end] =
                    MappabilityRun(last_run.start, base, simulated, correct, ectopic)
                continue
            end
        end
        push!(runs, MappabilityRun(base, base, simulated, correct, ectopic))
    end
    return runs
end

"""
Given one scaffold's runs and a span, expand the three counts to one value per
base of `[start_pos, stop_pos]`. Returns `(n_simulated, n_correct, n_ectopic)`
vectors, zero where no run covers a base.
"""
function expand_counts(
    runs::AbstractVector{MappabilityRun},
    start_pos::Integer,
    stop_pos::Integer,
)
    span = max(0, stop_pos - start_pos + 1)
    simulated, correct, ectopic =
        zeros(UInt16, span), zeros(UInt16, span), zeros(UInt32, span)
    for run in runs
        first_base, last_base = max(Int(run.start), start_pos), min(Int(run.stop), stop_pos)
        first_base > last_base && continue
        slots = (first_base-start_pos+1):(last_base-start_pos+1)
        simulated[slots] .= run.n_simulated
        correct[slots] .= run.n_correct
        ectopic[slots] .= run.n_ectopic
    end
    return (simulated, correct, ectopic)
end

#= Cache I/O =#

"""Name of the cache's manifest, holding its provenance and grid."""
const MANIFEST_FILE = "manifest.toml"

"""Name of the cache's scaffold table (scaffold name, file stem)."""
const SCAFFOLD_FILE = "scaffolds.tsv"

"""Identifies the cache layout; bump it when the layout changes."""
const CACHE_FORMAT = "BioinfoTools2 mappability runs v1"

"""Given a channel, the cache subdirectory holding it, e.g. `L150_F300_q2`."""
channel_dirname(channel) = "L$(channel.L)_F$(channel.F)_q$(channel.mapq_cutoff)"

"""
    write_mappability_arrow(filepath, runs; compress = :zstd)

Write one scaffold's runs to an Arrow file as five columns, the struct-of-arrays
layout of [`MappabilityRun`](@ref). `compress` is `:zstd`, `:lz4` or `nothing`;
only an uncompressed file can be memory-mapped on reading.
"""
function write_mappability_arrow(
    filepath::AbstractString,
    runs::AbstractVector{MappabilityRun};
    compress::Union{Nothing,Symbol} = :zstd,
)
    Arrow.write(filepath, StructArray(runs); compress = compress)
    return filepath
end

"""
    read_mappability_arrow(filepath; materialize = false)

Read an Arrow file written by [`write_mappability_arrow`](@ref) into a
`StructArray{MappabilityRun}` over the Arrow columns. An uncompressed file is
memory-mapped; `materialize = true` copies the columns into `Vector`s.
"""
function read_mappability_arrow(filepath::AbstractString; materialize::Bool = false)
    table = Arrow.Table(filepath)
    names = fieldnames(MappabilityRun)
    issubset(names, propertynames(table)) || throw(
        ArgumentError(
            "\"$filepath\" is not a mappability Arrow file (columns $(propertynames(table)))",
        ),
    )
    columns = map(
        name ->
            materialize ? collect(getproperty(table, name)) : getproperty(table, name),
        names,
    )
    return StructArray{MappabilityRun}(NamedTuple{names}(columns))
end

"""
Given a cache directory, read its `manifest.toml`. Returns the parsed `Dict`,
throwing when the directory is not a complete cache of this format.
"""
function mappability_manifest(dir::AbstractString)
    path = joinpath(dir, MANIFEST_FILE)
    isfile(path) ||
        throw(ArgumentError("no mappability cache at \"$dir\" (no $MANIFEST_FILE)"))
    manifest = TOML.parsefile(path)
    get(manifest, "format", "") == CACHE_FORMAT || throw(
        ArgumentError(
            "\"$dir\" is a cache of another format: $(get(manifest, "format", "?"))",
        ),
    )
    return manifest
end

"""
    mappability_channels(dir)

Given a cache directory, list the channels converted into it. Returns a
`Vector` of `(L, F, mapq_cutoff)` named tuples.
"""
mappability_channels(dir::AbstractString) = [
    MappabilityChannel((Int(L), Int(F), Int(cutoff))) for
    (L, F, cutoff) in mappability_manifest(dir)["channels"]
]

"""Given a cache directory, read its scaffold table. Returns `scaffold => stem` pairs."""
function scaffold_table(dir::AbstractString)
    pairs = Pair{String,String}[]
    for line in eachline(joinpath(dir, SCAFFOLD_FILE))
        tab = findfirst('\t', line)
        tab === nothing && continue
        push!(
            pairs,
            String(SubString(line, 1, tab - 1)) => String(SubString(line, tab + 1)),
        )
    end
    return pairs
end

"""
    read_mappability(dir; L, F, mapq_cutoff, materialize = false)

Given a cache written by [`convert_unimap`](@ref), open one channel. Returns a
[`MappabilityData`](@ref) whose scaffolds are memory-mapped when the cache is
uncompressed (decompressed into memory otherwise), copied when `materialize`.
"""
function read_mappability(
    dir::AbstractString;
    L::Integer,
    F::Integer,
    mapq_cutoff::Integer,
    materialize::Bool = false,
)
    manifest = mappability_manifest(dir)
    channel = MappabilityChannel((Int(L), Int(F), Int(mapq_cutoff)))
    channel in mappability_channels(dir) || throw(
        ArgumentError(
            "channel $channel is not in the cache; see mappability_channels(dir)",
        ),
    )
    channel_dir = joinpath(dir, channel_dirname(channel))
    scaffolds = Dict{String,StructArray{MappabilityRun}}(
        scaffold =>
            read_mappability_arrow(joinpath(channel_dir, "$stem.arrow"); materialize)
        for (scaffold, stem) in scaffold_table(dir)
    )
    return MappabilityData(scaffolds, channel, Float64(manifest["log_ratio_pseudocount"]))
end

"""
Given the scaffold names to be written, assign each a file stem. Returns
`scaffold => stem` pairs in sorted scaffold order (see `Methylation.sanitize_name`).
"""
function scaffold_stems(scaffolds)
    used = Set{String}()
    return [
        scaffold => _unique_stem!(used, sanitize_name(scaffold)) for
        scaffold in sort!(collect(scaffolds))
    ]
end

"""
Given a cache directory, its `scaffold => stem` pairs and the `manifest`, write
the scaffold table and then the manifest — last, so an interrupted conversion
leaves no manifest and [`read_mappability`](@ref) refuses the cache.
"""
function write_cache_index(dir::AbstractString, stems, manifest::AbstractDict)
    open(joinpath(dir, SCAFFOLD_FILE), "w") do io
        foreach(((scaffold, stem),) -> println(io, scaffold, '\t', stem), stems)
    end
    open(io -> TOML.print(io, manifest), joinpath(dir, MANIFEST_FILE), "w")
    return dir
end

"""
    convert_unimap(zarr_path, cache_dir, genome; channels, loci = nothing,
                   compress = :zstd, overwrite = false)

Given a UniMap `mappability.zarr` and the `Genome` it was built against, convert
`channels` — `(L, F, mapq_cutoff)` tuples, or `:all` — into a run cache at
`cache_dir`; a cache of the same store there gains the channels it lacks. **Needs
`using Zarr`.** `loci` (UniMap's `loci.tsv`) makes the genome check exact. Returns
`cache_dir`. **NOTE:** each call reads the whole store; see `docs/mappability.md` for sizes.
"""
function convert_unimap end

#= Masking and summaries =#

"""
    mask_calls(calls, data; below = 0.9, metric = :map_back, unmeasured = :drop)

Given methylation `calls`, drop those at bases whose `metric` is under `below`.
Calls on bases no run covers are dropped too, unless `unmeasured = :keep`.
Returns a new `MethylationData`. Each scaffold walks its calls and runs together.
"""
function mask_calls(
    calls::MethylationData,
    data::MappabilityData;
    below::Real = 0.9,
    metric::Symbol = :map_back,
    unmeasured::Symbol = :drop,
)
    unmeasured in (:drop, :keep) ||
        throw(ArgumentError("`unmeasured` must be :drop or :keep"))
    value_of = metric_function(metric, data.pseudocount)
    keep_unmeasured = unmeasured === :keep
    scaffold_names = sort!(collect(keys(calls.scaffolds)))
    tasks = map(scaffold_names) do scaffold
        Threads.@spawn begin
            scaffold_calls = calls.scaffolds[scaffold]
            runs = get(data.scaffolds, scaffold, nothing)
            keep =
                _keep_mask(scaffold_calls.pos, runs, value_of, below, keep_unmeasured)
            scaffold => scaffold_calls[keep]
        end
    end
    return MethylationData(Dict{String,StructArray{AggregatedCall}}(fetch.(tasks)))
end

"""
Given sorted call positions and one scaffold's runs (or `nothing`), decide per
call whether it stays. Returns a `BitVector`; one forward pass over both.
"""
function _keep_mask(
    positions,
    runs,
    value_of::V,
    below::Real,
    keep_unmeasured::Bool,
) where {V}
    keep = falses(length(positions))
    isnothing(runs) && return keep_unmeasured ? .!keep : keep
    entry = 1
    for (index, position) in enumerate(positions)
        while entry <= length(runs) && runs.stop[entry] < position
            entry += 1
        end
        covered = entry <= length(runs) && runs.start[entry] <= position
        keep[index] = covered ? !(value_of(runs[entry]) < below) : keep_unmeasured
    end
    return keep
end

"""
    binned_mappability(data; bin_size = 1000, metric = :map_back)

Given one channel, average `metric` over fixed bins of each scaffold, weighting
each base by its `n_simulated`. Returns `scaffold => Vector{Float32}` (bin `k`
covers bases `(k - 1) * bin_size + 1 : k * bin_size`), `NaN` where nothing was measured.
**NOTE:** scaffold lengths aren't known, so the bins stop at the last run.
"""
function binned_mappability(
    data::MappabilityData;
    bin_size::Integer = 1000,
    metric::Symbol = :map_back,
)
    bin_size >= 1 || throw(ArgumentError("`bin_size` must be positive"))
    value_of = metric_function(metric, data.pseudocount)
    binned = Dict{String,Vector{Float32}}()
    for (scaffold, runs) in data.scaffolds
        n_bins = isempty(runs) ? 0 : cld(Int(runs.stop[end]), bin_size)
        sums, weights = zeros(Float64, n_bins), zeros(Float64, n_bins)
        for run in runs
            value = value_of(run)
            isfinite(value) || continue
            for bin = cld(Int(run.start), bin_size):cld(Int(run.stop), bin_size)
                overlap =
                    min(Int(run.stop), bin * bin_size) -
                    max(Int(run.start), (bin - 1) * bin_size + 1) + 1
                sums[bin] += value * run.n_simulated * overlap
                weights[bin] += run.n_simulated * overlap
            end
        end
        binned[scaffold] = Float32.(
            map((total, weight) -> weight == 0 ? NaN : total / weight, sums, weights),
        )
    end
    return binned
end

#= Base.show overloads =#

function Base.show(io::IO, run::MappabilityRun)
    print(
        io,
        "MappabilityRun(",
        Int(run.start),
        ":",
        Int(run.stop),
        ", ",
        Int(run.n_correct),
        "/",
        Int(run.n_simulated),
        " correct, ",
        Int(run.n_ectopic),
        " ectopic)",
    )
end

function Base.show(io::IO, data::MappabilityData)
    (; L, F, mapq_cutoff) = data.channel
    print(
        io,
        "MappabilityData(L=$L F=$F MAPQ≥$mapq_cutoff, $(length(data)) scaffold",
        length(data) == 1 ? "" : "s",
        ", $(n_runs(data)) runs)",
    )
end

export MappabilityRun,
    MappabilityChannel,
    MappabilityData,
    append_runs!,
    binned_mappability,
    convert_unimap,
    ectopic_rate,
    expand_counts,
    log2_ratio,
    map_back,
    mappability_at,
    mappability_channels,
    mappability_manifest,
    mask_calls,
    metric_function,
    n_runs,
    per_base,
    read_mappability,
    read_mappability_arrow,
    runs_in_range,
    write_mappability_arrow

end
