"""
`convert_unimap`: UniMap `mappability.zarr` → BioinfoTools2's run cache.

An extension so that Zarr.jl (and its S3/GCS/HTTP backends) is loaded only by
those who convert. Reading the cache needs nothing from here.
"""
module BioinfoTools2ZarrExt

using IntervalTrees
using TOML
using Zarr

using BioinfoTools2.Reference: Genome, get_feature, get_strand, parse_strand
using BioinfoTools2.Data.Mappability:
    CACHE_FORMAT,
    MANIFEST_FILE,
    MappabilityChannel,
    MappabilityRun,
    append_runs!,
    channel_dirname,
    expand_counts,
    scaffold_stems,
    scaffold_table,
    write_cache_index,
    write_mappability_arrow
import BioinfoTools2.Data.Mappability: convert_unimap

"""Store rows read per call; caps a read of a Mb-scale locus at ~250 MB."""
const BLOCK_ROWS = 262_144

"""
Given the store's attributes, list its loci. Returns named tuples `(id, chrom,
start0, end0, strand, offset)`; locus `k`'s bases are rows `offset+1 : offset+end0-start0`.
"""
_store_loci(attrs) = [
    (
        id = String(locus["locus_id"]),
        chrom = String(locus["chrom"]),
        start0 = Int(locus["start0"]),
        end0 = Int(locus["end0"]),
        strand = only(String(locus["strand"])),
        offset = Int(locus["offset"]),
    ) for locus in attrs["loci"]
]

"""
Given UniMap's `loci.tsv` and the store's loci, read each locus's gene span.
Returns `locus_id => (gene_start0, gene_end0, strand)`, throwing when a store
locus is missing from the file — it belongs to another run.
"""
function _read_gene_spans(path::AbstractString, store_loci)
    lines = readlines(path)
    header = split(first(lines), '\t')
    column(name) = something(findfirst(==(name), header), 0)
    columns = column.(("locus_id", "gene_start0", "gene_end0", "strand"))
    any(iszero, columns) &&
        throw(ArgumentError("\"$path\" is not a UniMap loci.tsv (header: $header)"))
    spans =
        Dict{String,NamedTuple{(:gene_start0, :gene_end0, :strand),Tuple{Int,Int,Char}}}()
    for line in Iterators.drop(lines, 1)
        isempty(line) && continue
        fields = split(line, '\t')
        id, start0, end0, strand = getindex.(Ref(fields), columns)
        spans[id] = (
            gene_start0 = parse(Int, start0),
            gene_end0 = parse(Int, end0),
            strand = only(strand),
        )
    end
    absent = [locus.id for locus in store_loci if !haskey(spans, locus.id)]
    isempty(absent) || throw(
        ArgumentError(
            "$(length(absent)) store loci are not in \"$path\", e.g. $(first(absent))",
        ),
    )
    return spans
end

"""
Given one scaffold's gene tree and a locus, decide whether a `genome` gene
matches it: the exact gene span and strand from `loci.tsv` when `gene_spans` is
given, otherwise any same-strand gene inside the locus span.
"""
function _has_gene(tree, locus, gene_spans, minus_code)
    if isnothing(gene_spans)
        first_base, last_base, strand = locus.start0 + 1, locus.end0, locus.strand
        fits = interval -> interval.first >= first_base && interval.last <= last_base
    else
        span = gene_spans[locus.id]
        first_base, last_base, strand = span.gene_start0 + 1, span.gene_end0, span.strand
        fits = interval -> interval.first == first_base && interval.last == last_base
    end
    for interval in IntervalTrees.intersect(tree, (UInt32(first_base), UInt32(last_base)))
        fits(interval) &&
            (parse_strand(interval.value) == minus_code) == (strand == '-') &&
            return true
    end
    return false
end

"""
Given `genome` and the store's loci, check they belong together. Returns
`(n_loci, n_matched, method)`; throws when a locus scaffold is missing from
`genome` or no locus matches a gene, warns when only some do.
"""
function _match_genome(genome::Genome, store_loci, gene_spans)
    missing_scaffolds =
        sort!(unique(l.chrom for l in store_loci if !haskey(genome.scaffolds, l.chrom)))
    isempty(missing_scaffolds) || throw(
        ArgumentError(
            "$(length(missing_scaffolds)) locus scaffold(s) are not in `genome`, e.g. \"$(first(missing_scaffolds))\": wrong genome?",
        ),
    )
    genes = get_feature(genome, :gene)
    minus_code = get_strand('-')
    method =
        isnothing(gene_spans) ? "gene inside locus span, same strand" :
        "gene span and strand (loci.tsv)"
    unmatched = [
        locus.id for locus in store_loci if !(
            haskey(genes, locus.chrom) &&
            _has_gene(genes[locus.chrom], locus, gene_spans, minus_code)
        )
    ]
    n_matched = length(store_loci) - length(unmatched)
    n_matched == 0 && throw(
        ArgumentError(
            "none of the $(length(store_loci)) loci matches a gene of `genome` ($method): wrong genome or annotation?",
        ),
    )
    isempty(unmatched) ||
        @warn "$(length(unmatched)) of $(length(store_loci)) loci match no gene of `genome` ($method)" examples =
            first(unmatched, 5)
    return (n_loci = length(store_loci), n_matched = n_matched, method = method)
end

"""
Given the requested `channels` (`:all`, or `(L, F, mapq_cutoff)` tuples) and the
store's grid, resolve them to `(iL, iF, ic)` array indices, in request order.
"""
function _channel_indices(channels, read_lengths, fragment_lengths, cutoffs)
    channels === :all && return [
        (iL, iF, ic) for iL in eachindex(read_lengths) for
        iF in eachindex(fragment_lengths) for ic in eachindex(cutoffs)
    ]
    indices = map(collect(channels)) do (L, F, cutoff)
        found = (
            findfirst(==(L), read_lengths),
            findfirst(==(F), fragment_lengths),
            findfirst(==(cutoff), cutoffs),
        )
        any(isnothing, found) &&
            throw(ArgumentError("channel ($L, $F, $cutoff) is not on the store's grid"))
        found
    end
    return unique(indices)
end

"""
Given the store's three count arrays and a range of store rows, read them for
every channel at once — each Zarr chunk holds all channels anyway.
"""
_read_rows(arrays, rows) = (
    n_simulated = arrays.n_simulated[rows, :, :],
    n_correct = arrays.n_correct[rows, :, :, :],
    n_ectopic = arrays.n_ectopic[rows, :, :, :],
)

"""Given a locus and a 0-based genomic span `[start0, end0)` inside it, its store rows."""
_locus_rows(locus, start0::Int, end0::Int) =
    (locus.offset+start0-locus.start0+1):(locus.offset+end0-locus.start0)

"""
Given per-channel runs and a block of counts starting at 1-based `first_base`,
extend every channel's runs with it, one task per channel.
"""
function _encode_block!(runs, block, indices, first_base::Int)
    Threads.@threads for slot in eachindex(indices)
        iL, iF, ic = indices[slot]
        append_runs!(
            runs[slot],
            first_base,
            view(block.n_simulated, :, iF, iL),
            view(block.n_correct, :, iF, ic, iL),
            view(block.n_ectopic, :, iF, ic, iL),
        )
    end
    return runs
end

"""
Given the runs already written and a locus overlapping them on `[start0, end0)`,
check the locus holds the same counts there (UniMap cuts both from one region
track, so they must). Throws on any difference.
"""
function _check_overlap(runs, arrays, indices, locus, start0::Int, end0::Int)
    block = _read_rows(arrays, _locus_rows(locus, start0, end0))
    probe = MappabilityRun(0, start0 + 1, 0, 0, 0)
    agrees = map(eachindex(indices)) do slot
        iL, iF, ic = indices[slot]
        tail = view(
            runs[slot],
            searchsortedfirst(runs[slot], probe; by = run->run.stop):length(runs[slot]),
        )
        expected = expand_counts(tail, start0 + 1, end0)
        expected == (
            block.n_simulated[:, iF, iL],
            block.n_correct[:, iF, ic, iL],
            block.n_ectopic[:, iF, ic, iL],
        )
    end
    all(agrees) || error(
        "locus $(locus.id) disagrees with an overlapping locus on $(locus.chrom):$(start0 + 1)-$end0",
    )
    return nothing
end

"""
Given the count arrays, one scaffold's loci and the channels, run-encode the
scaffold. Returns one `Vector{MappabilityRun}` per channel. Loci go in genomic
order; a base shared with an earlier locus is checked, not encoded twice.
"""
function _convert_scaffold(arrays, scaffold_loci, indices)
    runs = [MappabilityRun[] for _ in indices]
    encoded_end = 0   # 0-based exclusive end of what is already encoded
    for locus in sort(scaffold_loci; by = locus -> (locus.start0, locus.end0))
        overlap_end = min(locus.end0, encoded_end)
        overlap_end > locus.start0 &&
            _check_overlap(runs, arrays, indices, locus, locus.start0, overlap_end)
        for block_start = max(locus.start0, encoded_end):BLOCK_ROWS:(locus.end0-1)
            block_end = min(block_start + BLOCK_ROWS, locus.end0)
            block = _read_rows(arrays, _locus_rows(locus, block_start, block_end))
            _encode_block!(runs, block, indices, block_start + 1)
        end
        encoded_end = max(encoded_end, locus.end0)
    end
    return runs
end

"""
Given a UniMap store's path, open it and read what conversion needs. Returns
`(store, attrs, grid, loci)`, `grid` holding `read_lengths`, `fragment_lengths`
and `mapq_cutoffs`; throws when the store has no count arrays.
"""
function _open_store(zarr_path::AbstractString)
    store = zopen(zarr_path; consolidated = isfile(joinpath(zarr_path, ".zmetadata")))
    all(name -> haskey(store.arrays, name), ("n_simulated", "n_correct", "n_ectopic")) ||
        throw(ArgumentError("\"$zarr_path\" has no UniMap count arrays"))
    attrs = store.attrs
    grid = (
        read_lengths = Int.(attrs["read_lengths"]),
        fragment_lengths = Int.(attrs["fragment_lengths"]),
        mapq_cutoffs = Int.(attrs["mapq_cutoffs"]),
    )
    return (store, attrs, grid, _store_loci(attrs))
end

"""
Given the output directory and the store, decide what to build on. Returns the
existing cache's manifest when channels are being added to a cache of this same
store, or `nothing` once the directory is ready for a fresh cache (cleared first
when `overwrite`). Refuses to touch a directory that is not a mappability cache.
"""
function _existing_cache(dir::AbstractString, attrs, grid, overwrite::Bool)
    (isdir(dir) && !isempty(readdir(dir))) || (mkpath(dir); return nothing)
    manifest_path = joinpath(dir, MANIFEST_FILE)
    manifest = isfile(manifest_path) ? TOML.parsefile(manifest_path) : Dict{String,Any}()
    get(manifest, "format", "") == CACHE_FORMAT || throw(
        ArgumentError(
            "\"$dir\" is not empty and not a mappability cache; refusing to write there",
        ),
    )
    if overwrite
        rm(dir; recursive = true)
        mkpath(dir)
        return nothing
    end
    same_store =
        manifest["genome_checksum"] == string(get(attrs, "genome_checksum", "")) &&
        all(key -> manifest[string(key)] == grid[key], keys(grid)) &&
        manifest["n_loci"] == length(attrs["loci"])
    same_store || throw(
        ArgumentError(
            "\"$dir\" caches a different store; pass overwrite = true to replace it",
        ),
    )
    return manifest
end

"""
Given the store's arrays, its loci by scaffold, the scaffold file stems and the
channels, convert every scaffold and write one Arrow file per scaffold per
channel. Returns each channel's total run count.
"""
function _write_channels(
    cache_dir,
    arrays,
    by_scaffold,
    stems,
    indices,
    channel_list,
    compress,
)
    foreach(channel -> mkpath(joinpath(cache_dir, channel_dirname(channel))), channel_list)
    run_counts = zeros(Int, length(indices))
    for (scaffold, stem) in stems
        runs = _convert_scaffold(arrays, by_scaffold[scaffold], indices)
        for (slot, channel) in enumerate(channel_list)
            path = joinpath(cache_dir, channel_dirname(channel), "$stem.arrow")
            write_mappability_arrow(path, runs[slot]; compress)
            run_counts[slot] += length(runs[slot])
        end
    end
    return run_counts
end

"""
Given the store, a fresh cache's manifest. Returns the `Dict`, with no channels
yet: [`convert_unimap`](@ref) appends them as they are converted.
"""
_new_manifest(zarr_path, attrs, grid, store_loci, genome_match, compress) =
    Dict{String,Any}(
        "format" => CACHE_FORMAT,
        "source" => abspath(zarr_path),
        "genome_checksum" => string(get(attrs, "genome_checksum", "")),
        "read_lengths" => grid.read_lengths,
        "fragment_lengths" => grid.fragment_lengths,
        "mapq_cutoffs" => grid.mapq_cutoffs,
        "n_loci" => length(store_loci),
        "log_ratio_pseudocount" => Float64(get(attrs, "log_ratio_pseudocount", 0.5)),
        "compress" => isnothing(compress) ? "none" : string(compress),
        "genome_match" =>
            Dict(string(key) => value for (key, value) in pairs(genome_match)),
        "channels" => Vector{Int}[],
        "n_runs" => Int[],
    )

function convert_unimap(
    zarr_path::AbstractString,
    cache_dir::AbstractString,
    genome::Genome;
    channels = nothing,
    loci::Union{Nothing,AbstractString} = nothing,
    compress::Union{Nothing,Symbol} = :zstd,
    overwrite::Bool = false,
)
    isnothing(channels) && throw(
        ArgumentError(
            "name the channels to convert, e.g. channels = [(150, 300, 2)], or channels = :all " *
            "(~0.1-1.3 GB each for a 1 Gb panel; see docs/mappability.md)",
        ),
    )
    store, attrs, grid, store_loci = _open_store(zarr_path)
    gene_spans = isnothing(loci) ? nothing : _read_gene_spans(loci, store_loci)
    genome_match = _match_genome(genome, store_loci, gene_spans)
    indices = _channel_indices(
        channels,
        grid.read_lengths,
        grid.fragment_lengths,
        grid.mapq_cutoffs,
    )

    manifest = _existing_cache(cache_dir, attrs, grid, overwrite)
    if isnothing(manifest)
        manifest = _new_manifest(zarr_path, attrs, grid, store_loci, genome_match, compress)
        by_chrom = unique(locus.chrom for locus in store_loci)
        stems = scaffold_stems(by_chrom)
    else
        # Adding to a cache: same file names, same compression, only new channels.
        stems = scaffold_table(cache_dir)
        compress = manifest["compress"] == "none" ? nothing : Symbol(manifest["compress"])
        filter!(
            ((iL, iF, ic),) -> !(
                [grid.read_lengths[iL], grid.fragment_lengths[iF], grid.mapq_cutoffs[ic]] in manifest["channels"]
            ),
            indices,
        )
        isempty(indices) && return cache_dir
    end
    channel_list = [
        MappabilityChannel((
            grid.read_lengths[iL],
            grid.fragment_lengths[iF],
            grid.mapq_cutoffs[ic],
        )) for (iL, iF, ic) in indices
    ]

    by_scaffold = Dict{String,Vector{eltype(store_loci)}}()
    foreach(
        locus -> push!(get!(by_scaffold, locus.chrom, eltype(store_loci)[]), locus),
        store_loci,
    )
    arrays = (
        n_simulated = store["n_simulated"],
        n_correct = store["n_correct"],
        n_ectopic = store["n_ectopic"],
    )
    run_counts = _write_channels(
        cache_dir,
        arrays,
        by_scaffold,
        stems,
        indices,
        channel_list,
        compress,
    )

    append!(
        manifest["channels"],
        [[channel.L, channel.F, channel.mapq_cutoff] for channel in channel_list],
    )
    append!(manifest["n_runs"], run_counts)
    manifest["updated"] = Libc.strftime("%Y-%m-%dT%H:%M:%S", time())
    write_cache_index(cache_dir, stems, manifest)
    return cache_dir
end

end
