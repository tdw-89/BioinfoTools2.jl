module Reference

using BioGenerics
using CodecZlib
using GFF3
using IntervalTrees

using ..BitCodes
using ..SOTerms

"""How many parsed records `add_features!` batches onto its channel at a time."""
const RECORD_BUFFER = 1_000

"""Interval tree of 1-based closed `[start, end]` spans carrying a 64-bit code."""
const IntervalTreeM64 = IntervalTree{UInt32,IntervalValue{UInt32,UInt64}}

"""
One sequence of a genome:
- `name`: the name of the scaffold
- `features`: an interval tree of feature spans, each carrying a 64-bit
  metadata code (see [`Genome`](@ref)).
"""
struct Scaffold
    name::String
    features::IntervalTreeM64
end

"""
A bare `(start_pos, end_pos, code)` interval, for callers that need no scaffold
or vocabulary context (e.g. `Homologs.Paralogs.ParalogGroup`'s per-scaffold
`StructArray`s).
"""
struct IntervalSimple
    start_pos::UInt32
    end_pos::UInt32
    code::UInt64
end

"""One GFF3 record reduced to what [`add_features!`](@ref) commits to a `Genome`."""
struct ParseResult
    scaffold_id::String
    start_pos::UInt32
    end_pos::UInt32
    code::UInt64
    id::String
    source::String
    biotype::String
end

"""
The top level of an in-memory genome representation: the `scaffolds` keyed by
name, a `vocab`/`vocab_lookup` string intern pool for feature metadata, and the
`meta_offsets`/`meta_blob` byte store the intervals index into.

Each feature's metadata is split between the 64-bit code on its interval and
that blob. Bits 0-31 of the code are a 1-based index into `meta_offsets`, giving
the feature's start offset in `meta_blob`; the remaining bits hold the strand
and SO term (see [`pack_metadata`](@ref)). A blob entry is a run of `UInt32`
vocabulary tokens — ID, source and biotype — resolved through `vocab`.
"""
mutable struct Genome
    scaffolds::Dict{String,Scaffold}

    # String intern pool
    vocab::Vector{String}
    vocab_lookup::Dict{String,UInt32}

    # Metadata store
    meta_offsets::Vector{UInt32}
    meta_blob::Vector{UInt8}
end

"""
A single annotated feature looked up from a `Genome`:
- `id`: the ID the feature was searched by (first metadata field)
- `feature_type`: the SO term label (e.g. `:gene`)
- `chromosome`: the scaffold name the feature lives on
- `start_pos` / `end_pos`: 1-based closed interval bounds
- `metadata`: the parsed metadata stored for the feature (ID, source, biotype)
- `code`: the raw 64-bit metadata code (strand, SO term, metadata index)
"""
struct FeatureRecord
    id::String
    feature_type::Symbol
    chromosome::String
    start_pos::UInt32
    end_pos::UInt32
    metadata::Vector{String}
    code::UInt64
end

"""A named organism and its [`Genome`](@ref)."""
mutable struct Species
    name::String
    taxon_id::String
    genome::Genome
end

"""Create a `Species` with an empty genome."""
function Species(name::String; taxon_id::String = "")
    genome =
        Genome(Dict{String,Scaffold}(), String[], Dict{String,UInt32}(), UInt32[], UInt8[])
    return Species(name, taxon_id, genome)
end

"""Get the 16-bit code for a given SO term label, or `nothing` if unknown."""
function convert_so_term(label::AbstractString)
    result = SO_TERMS[label]
    return isnothing(result) ? nothing : result[1]
end

"""
Strip common feature-type prefixes (for example `gene:` or `transcript:`)
from IDs when present. Prefixes can be stacked (`gene:transcript:...`).
"""
function sanitize_id(id::AbstractString)
    cleaned = String(id)
    # Every prefix ends in ':', so an ID without one skips the regex.
    while occursin(':', cleaned)
        matched = match(r"^(gene|transcript|mrna|rna|cds|protein):(.+)$"i, cleaned)
        matched === nothing && return cleaned
        cleaned = String(matched.captures[2])
    end
    return cleaned
end

"""The GFF3 attribute a feature's biotype is read from (NCBI's convention)."""
const BIOTYPE_KEY = "gene_biotype"

"""
Given a string and a byte range of `data`, compare them without allocating.
Returns `true` when the bytes are exactly the string's.
"""
@inline function _bytes_equal(text::String, data::Vector{UInt8}, range::UnitRange{Int})
    sizeof(text) == length(range) || return false
    isempty(range) && return true
    return GC.@preserve text data ccall(
        :memcmp,
        Cint,
        (Ptr{UInt8}, Ptr{UInt8}, Csize_t),
        pointer(text),
        pointer(data, first(range)),
        length(range),
    ) == 0
end

"""
Single value of GFF3 attribute `key` on `record`, or `"NA"` when absent or
multi-valued. Only that attribute is decoded; a repeated key reads its last copy.
"""
function attribute_or_na(record::GFF3.Record, key::String)
    found = 0
    for (index, key_range) in enumerate(record.attribute_keys)
        _bytes_equal(key, record.data, key_range) && (found = index)
    end
    found == 0 && return "NA"
    values = GFF3.attrvals(record, found)
    return length(values) == 1 ? only(values) : "NA"
end

"""
Each distinct type column one load has met, as raw bytes, beside the SO short
code it loads under: `0` when the type is not a SOFA term or not in `wanted`
(`nothing` wants every term). Saves resolving `SO_TERMS` once per record.
"""
struct FeatureTypeCache
    labels::Vector{String}
    codes::Vector{UInt16}
    wanted::Union{Nothing,Set{UInt16}}
end

"""Create an empty `FeatureTypeCache` keeping the SO codes in `wanted` (`nothing`: all)."""
FeatureTypeCache(wanted::Union{Nothing,Set{UInt16}} = nothing) =
    FeatureTypeCache(String[], UInt16[], wanted)

"""
Given a record, resolve its type column through `cache`, adding the type on a
miss. Returns the SO short code to load it under, `0` to skip it.
"""
function _so_code!(cache::FeatureTypeCache, record::GFF3.Record)
    for (slot, label) in enumerate(cache.labels)
        _bytes_equal(label, record.data, record.type_) && return cache.codes[slot]
    end
    code = something(convert_so_term(GFF3.featuretype(record)), 0x0000)
    wanted = cache.wanted
    isnothing(wanted) || code in wanted || (code = 0x0000)
    push!(cache.labels, String(record.data[record.type_]))
    push!(cache.codes, code)
    return code
end

"""
The last value decoded from one GFF3 column, under the raw bytes it came from,
so a run of records on one scaffold (or from one source) shares one `String`.
"""
mutable struct ColumnCache
    raw::String
    value::String
end

"""Create an empty `ColumnCache`: no column holds a tab, so its first lookup misses."""
ColumnCache() = ColumnCache("\t", "")

"""
Given a record and its column's byte `range`, return the cached value when the
bytes match, else `read_value(record)`, caching it.
"""
@inline function _column!(read_value, cache::ColumnCache, record::GFF3.Record, range)
    _bytes_equal(cache.raw, record.data, range) && return cache.value
    cache.value = read_value(record)
    cache.raw = String(record.data[range])
    return cache.value
end

"""The source column of `record`, `"NA"` when it is `.`."""
_source_or_na(record::GFF3.Record) = GFF3.hassource(record) ? GFF3.source(record) : "NA"

"""
Everything one load decides up front — which SO codes and biotypes to keep
(`nothing` keeps all), whether to sanitize IDs — plus its column caches.
"""
struct RecordParser
    sanitize_ids::Bool
    biotypes::Union{Nothing,Set{String}}
    types::FeatureTypeCache
    scaffolds::ColumnCache
    sources::ColumnCache
end

"""
Create a `RecordParser`, taking `feature_types` and `biotypes` as
`_wanted_types`/`_wanted_biotypes` return them.
"""
RecordParser(;
    sanitize_ids::Bool = true,
    feature_types::Union{Nothing,Set{UInt16}} = nothing,
    biotypes::Union{Nothing,Set{String}} = nothing,
) = RecordParser(
    sanitize_ids,
    biotypes,
    FeatureTypeCache(feature_types),
    ColumnCache(),
    ColumnCache(),
)

"""
Given a record, reduce it to a [`ParseResult`](@ref) carrying `meta_index`.
Returns `nothing` for a type or biotype `parser` does not keep; the type is
checked first, so a skipped type costs no attribute decoding.
"""
function _parse!(parser::RecordParser, record::GFF3.Record, meta_index::UInt32)
    so_term = _so_code!(parser.types, record)
    so_term == 0x0000 && return nothing
    biotype = attribute_or_na(record, BIOTYPE_KEY)
    wanted = parser.biotypes
    isnothing(wanted) || biotype in wanted || return nothing

    feature_id = attribute_or_na(record, "ID")
    parser.sanitize_ids && (feature_id = sanitize_id(feature_id))
    return ParseResult(
        _column!(GFF3.seqname, parser.scaffolds, record, record.seqid),
        UInt32(GFF3.seqstart(record)),
        UInt32(GFF3.seqend(record)),
        pack_metadata(meta_index, strand_code(GFF3.strand(record)), so_term),
        feature_id,
        _column!(_source_or_na, parser.sources, record, record.source),
        biotype,
    )
end

"""
Reduce a GFF3 record to a [`ParseResult`](@ref) carrying metadata index
`meta_index`, or `nothing` when its feature type is not a known SO term.
Feature IDs are stripped of type prefixes unless `sanitize_ids` is `false`.
"""
parse_record(record::GFF3.Record, meta_index::UInt32; sanitize_ids::Bool = true) =
    _parse!(RecordParser(; sanitize_ids), record, meta_index)

#= Metadata store =#

"""Read the `UInt32` vocabulary token stored at 1-based byte `offset` of `blob`."""
@inline _read_token(blob::Vector{UInt8}, offset::Integer) =
    GC.@preserve blob unsafe_load(Ptr{UInt32}(pointer(blob, offset)))

"""
Given a vocabulary token, append its 4 bytes to `blob` (native endian, matching
[`_read_token`](@ref)). Returns `blob`.
"""
@inline function _append_token!(blob::Vector{UInt8}, token::UInt32)
    offset = length(blob) + 1
    resize!(blob, offset + 3)
    GC.@preserve blob unsafe_store!(Ptr{UInt32}(pointer(blob, offset)), token)
    return blob
end

"""
Given a metadata index, locate its entry in the blob. Returns the byte range, or
`nothing` when the index is out of range or the entry is empty.
"""
@inline function _meta_range(genome::Genome, meta_index::UInt32)
    index = Int(meta_index)
    (index < 1 || index + 1 > length(genome.meta_offsets)) && return nothing
    start_byte = Int(genome.meta_offsets[index])
    end_byte = Int(genome.meta_offsets[index+1]) - 1
    return start_byte > end_byte ? nothing : (start_byte:end_byte)
end

"""
Resolve the interned metadata strings (ID, source, biotype) stored for
`meta_index`. Returns an empty vector when the index carries no metadata.
"""
function get_metadata(genome::Genome, meta_index::UInt32)
    bytes = _meta_range(genome, meta_index)
    bytes === nothing && return String[]
    blob = genome.meta_blob
    return [genome.vocab[_read_token(blob, offset)] for offset = first(bytes):4:last(bytes)]
end

"""
Return only the interned ID (the first metadata token) for `meta_index`, or
`nothing` when the feature carries no metadata.

Reads just the leading 4 bytes of the entry rather than allocating the full
metadata vector; prefer this over [`get_metadata`](@ref) when only the ID is
needed.
"""
function get_metadata_id(genome::Genome, meta_index::UInt32)
    bytes = _meta_range(genome, meta_index)
    bytes === nothing && return nothing
    return genome.vocab[_read_token(genome.meta_blob, first(bytes))]
end

"""
Given a feature interval, resolve the interned ID stored for it. Returns that
ID, or `nothing` when the feature carries no metadata.
"""
feature_id(genome::Genome, interval::IntervalValue{UInt32,UInt64}) =
    get_metadata_id(genome, parse_index(interval.value))

"""Metadata for every feature of `features`, in the tree's own order."""
function get_metadata(genome::Genome, features::IntervalTreeM64)
    return [get_metadata(genome, parse_index(interval.value)) for interval in features]
end

get_metadata(genome::Genome, scaffold::Scaffold) = get_metadata(genome, scaffold.features)

"""Metadata for every feature of every scaffold, keyed by scaffold name."""
function get_metadata(genome::Genome)
    return Dict{String,Vector{Vector{String}}}(
        name => get_metadata(genome, scaffold) for (name, scaffold) in genome.scaffolds
    )
end

#= Feature lookup =#

"""
Interval tree of the `feature`-type features of `scaffold`, or `nothing` when
`feature` is not a known SO term.
"""
function get_feature(scaffold::Scaffold, feature::Symbol)::Union{Nothing,IntervalTreeM64}
    result = SO_TERMS[feature]
    isnothing(result) && return nothing
    so_term = result[1]

    tree = IntervalTreeM64()
    for interval in scaffold.features
        parse_so_term(interval.value) == so_term && push!(tree, interval)
    end
    return tree
end

get_feature(scaffold::Scaffold, feature::AbstractString)::Union{Nothing,IntervalTreeM64} =
    get_feature(scaffold, Symbol(feature))

"""
The `feature`-type features of every scaffold, keyed by scaffold name. Returns
an empty `Dict` when `feature` is not a known SO term.
"""
function get_feature(genome::Genome, feature::Symbol)
    isnothing(SO_TERMS[feature]) && return Dict{String,IntervalTreeM64}()
    return Dict{String,IntervalTreeM64}(
        name => get_feature(scaffold, feature) for (name, scaffold) in genome.scaffolds
    )
end

get_feature(genome::Genome, feature::AbstractString) = get_feature(genome, Symbol(feature))

"""Return the unique SO term labels used by features stored in `genome`."""
function get_so_terms(genome::Genome)
    terms = Set{Symbol}()
    for scaffold in values(genome.scaffolds), interval in scaffold.features
        term = SO_TERMS[parse_so_term(interval.value)]
        isnothing(term) || push!(terms, term[2])
    end
    return sort!(collect(terms))
end

"""
Given a genome, lazily iterate `(scaffold_name, interval)` over every feature in
it. Returns that generator, the shared walk behind the `getindex` methods below.
"""
_features(genome::Genome) = (
    (name, interval) for (name, scaffold) in genome.scaffolds for
    interval in scaffold.features
)

"""Build a [`FeatureRecord`](@ref) from a scaffold name and one of its intervals."""
function feature_record(
    genome::Genome,
    chromosome::String,
    interval::IntervalValue{UInt32,UInt64},
)
    code = interval.value
    so_result = SO_TERMS[parse_so_term(code)]
    metadata = get_metadata(genome, parse_index(code))
    return FeatureRecord(
        isempty(metadata) ? "" : metadata[1],
        isnothing(so_result) ? Symbol("") : so_result[2],
        chromosome,
        interval.first,
        interval.last,
        metadata,
        code,
    )
end

"""
Look up a single feature by its ID (the first metadata field), returning the
first match as a [`FeatureRecord`](@ref) or `nothing`.

This scans the genome's features (O(features) per lookup) rather than keeping an
ID index; use the `Vector` method to resolve many IDs in one walk.
"""
function Base.getindex(genome::Genome, id::AbstractString)
    for (name, interval) in _features(genome)
        if feature_id(genome, interval) == id
            return feature_record(genome, name, interval)
        end
    end
    return nothing
end

"""
Look up every feature whose ID (first metadata field) is contained in `ids`,
returning a `Vector{FeatureRecord}` resolved in one O(features) walk. IDs with
no matching feature are skipped.
"""
function Base.getindex(genome::Genome, ids::AbstractVector{<:AbstractString})
    wanted = Set{String}(ids)
    records = FeatureRecord[]
    for (name, interval) in _features(genome)
        found = feature_id(genome, interval)
        if !isnothing(found) && found in wanted
            push!(records, feature_record(genome, name, interval))
        end
    end
    return records
end

"""
Look up the feature carrying the 32-bit metadata index `meta_index`, returning
it as a [`FeatureRecord`](@ref), or `nothing` when no feature uses that index.

Scans the genome's intervals (O(features) per lookup); use the `Vector` method
to resolve many indices in one walk.
"""
function Base.getindex(genome::Genome, meta_index::UInt32)
    for (name, interval) in _features(genome)
        if parse_index(interval.value) == meta_index
            return feature_record(genome, name, interval)
        end
    end
    return nothing
end

"""
Look up every feature whose 32-bit metadata index appears in `meta_indices`,
returning a `Dict` mapping each found index to its [`FeatureRecord`](@ref).
Indices matching no feature are absent from the result.

One O(features) walk resolves the whole set and stops early once every requested
index has been found, so callers holding many indices must use this rather than
the scalar method, which would be O(indices × features). As in the scalar
method, the first interval carrying an index wins.
"""
function Base.getindex(genome::Genome, meta_indices::AbstractVector{UInt32})
    records = Dict{UInt32,FeatureRecord}()
    wanted = Set{UInt32}(meta_indices)
    isempty(wanted) && return records

    for (name, interval) in _features(genome)
        meta_index = parse_index(interval.value)
        if meta_index in wanted
            records[meta_index] = feature_record(genome, name, interval)
            delete!(wanted, meta_index)
            isempty(wanted) && return records
        end
    end
    return records
end

#= Loading =#

"""
Given a metadata string, intern it into the genome's vocab. Returns its 1-based
`UInt32` token.
"""
function _intern_string!(genome::Genome, s::String)
    get!(genome.vocab_lookup, s) do
        token = UInt32(length(genome.vocab) + 1)
        push!(genome.vocab, s)
        token
    end
end

"""
Given a genome, the metadata index its next loaded feature takes. Returns one
past the features already stored (`meta_offsets` ends in a sentinel once loaded).
"""
_next_meta_index(genome::Genome) = UInt32(max(length(genome.meta_offsets), 1))

"""
Given a channel of parsed record batches, drain it into `genome` — intervals and
metadata blob alike. Returns `nothing`; runs on its own task, so that committing
overlaps [`add_features!`](@ref)'s parsing.
"""
function _build_genome!(channel::Channel{Vector{ParseResult}}, genome::Genome)
    # A previous load's closing sentinel is exactly this load's first offset.
    isempty(genome.meta_offsets) || pop!(genome.meta_offsets)
    scaffold = nothing
    for batch in channel, result in batch
        if isnothing(scaffold) || scaffold.name != result.scaffold_id
            scaffold = get!(genome.scaffolds, result.scaffold_id) do
                Scaffold(result.scaffold_id, IntervalTreeM64())
            end
        end
        push!(
            scaffold.features,
            IntervalValue(result.start_pos, result.end_pos, result.code),
        )

        push!(genome.meta_offsets, UInt32(length(genome.meta_blob) + 1))
        for field in (result.id, result.source, result.biotype)
            _append_token!(genome.meta_blob, _intern_string!(genome, field))
        end
    end

    # Sentinel offset so the last feature's end byte can be computed.
    push!(genome.meta_offsets, UInt32(length(genome.meta_blob) + 1))
    return nothing
end

"""
Given an open file handle and the path it came from, wrap it for reading.
Returns a decompressing stream when the name ends in `.gz`, and `handle`
unchanged otherwise.
"""
gzip_stream(handle::IO, path::AbstractString) =
    endswith(path, ".gz") ? GzipDecompressorStream(handle) : handle

"""
Given an open GFF3 file handle and its path, wrap it in a reader that keeps
directives and comments. Returns that `GFF3.Reader`; callers skip non-feature
records themselves.

**NOTE:** GFF3.jl consumes a skipped line silently, then throws `EOFError` when
only such lines (NCBI's closing `###`) follow the last record.
"""
_gff_reader(handle::IO, path::AbstractString) =
    GFF3.Reader(gzip_stream(handle, path); skip_directives = false, skip_comments = false)

"""
Given an open GFF3 reader, parse its records onto `channel` in batches of
[`RECORD_BUFFER`](@ref), numbering them from `meta_index`. Returns `nothing`;
directives, comments and records `parser` skips consume no metadata index.
"""
function _batch_records!(
    channel::Channel{Vector{ParseResult}},
    reader::GFF3.Reader,
    parser::RecordParser,
    meta_index::UInt32,
)
    record = GFF3.Record()
    buffer = sizehint!(ParseResult[], RECORD_BUFFER)

    while !eof(reader)
        # No need to `empty!` record; GFF3.jl does that inside `read!`
        read!(reader, record)
        (BioGenerics.isfilled(record) && GFF3.isfeature(record)) || continue

        result = _parse!(parser, record, meta_index)
        isnothing(result) && continue

        push!(buffer, result)
        meta_index += UInt32(1)
        if length(buffer) == RECORD_BUFFER
            put!(channel, buffer)
            buffer = sizehint!(ParseResult[], RECORD_BUFFER)
        end
    end
    # Flush any remaining records that didn't fill a complete batch
    isempty(buffer) || put!(channel, buffer)
    return nothing
end

"""
Given SO term labels (one, several, or `nothing`), resolve their short codes.
Returns a `Set{UInt16}`, or `nothing` to keep every type; throws an
`ArgumentError` on a label that is not a SOFA term.
"""
_wanted_types(::Nothing) = nothing
_wanted_types(label::Union{Symbol,AbstractString}) = _wanted_types((label,))
function _wanted_types(labels)
    codes = Set{UInt16}()
    for label in labels
        found = SO_TERMS[Symbol(label)]
        isnothing(found) && throw(ArgumentError("`$label` is not a SOFA feature type"))
        push!(codes, found[1])
    end
    return codes
end

"""
Given biotypes (one, several, or `nothing`), collect them. Returns a
`Set{String}`, or `nothing` to keep every biotype.
"""
_wanted_biotypes(::Nothing) = nothing
_wanted_biotypes(biotype::Union{Symbol,AbstractString}) = Set([String(biotype)])
_wanted_biotypes(biotypes) = Set{String}(String(biotype) for biotype in biotypes)

"""
    add_features!(gff_path, genome; sanitize_ids, feature_types, biotypes)

Given a GFF3 file (optionally gzipped), add its SOFA-typed features to `genome`.
Returns `nothing`. `feature_types` (SO labels, e.g. `:gene`) and `biotypes`
(`gene_biotype` values, e.g. `"protein_coding"`) each keep only matching records;
`nothing` keeps all. IDs lose type prefixes unless `sanitize_ids` is `false`.

**NOTE:** a record without `gene_biotype` (in NCBI files every `mRNA`, `exon` and
`CDS`) has biotype `"NA"`, so `biotypes` drops it unless `"NA"` is listed.
Types SOFA lacks (NCBI's `lncRNA`, `V_gene_segment`) never load, and naming one throws.
"""
function add_features!(
    gff_path::String,
    genome::Genome;
    sanitize_ids::Bool = true,
    feature_types = nothing,
    biotypes = nothing,
)
    parser = RecordParser(;
        sanitize_ids,
        feature_types = _wanted_types(feature_types),
        biotypes = _wanted_biotypes(biotypes),
    )
    # Read before the builder spawns: it pops the previous load's sentinel.
    first_index = _next_meta_index(genome)
    # Unbounded channel so the parser never blocks waiting for the builder
    channel = Channel{Vector{ParseResult}}(Inf)
    builder_task = Threads.@spawn _build_genome!(channel, genome)

    open(gff_path) do handle
        reader = _gff_reader(handle, gff_path)
        try
            _batch_records!(channel, reader, parser, first_index)
        finally
            close(reader)
            close(channel)
        end
    end

    wait(builder_task)
    return nothing
end

add_features!(gff_path::String, species::Species; options...) =
    add_features!(gff_path, species.genome; options...)

"""
Given a GFF3 file, read the links [`protein_gene_map`](@ref) climbs. Returns
`(genes, parents, cds_parents)`: the gene IDs, every other feature's `ID =>
Parent` (exons, which parent nothing, aside) and each CDS's `protein_id => Parent`.
"""
function _protein_links(gff_path::String, sanitize_ids::Bool)
    gene_code, cds_code, exon_code = (SO_TERMS[term][1] for term in (:gene, :CDS, :exon))
    types = FeatureTypeCache()
    clean(id) = sanitize_ids ? sanitize_id(id) : id
    genes, parents, cds_parents =
        Set{String}(), Dict{String,String}(), Dict{String,String}()

    open(gff_path) do handle
        reader = _gff_reader(handle, gff_path)
        record = GFF3.Record()
        try
            while !eof(reader)
                read!(reader, record)
                (BioGenerics.isfilled(record) && GFF3.isfeature(record)) || continue
                code = _so_code!(types, record)
                code == exon_code && continue
                id = attribute_or_na(record, "ID")
                if code == gene_code
                    id == "NA" || push!(genes, clean(id))
                    continue
                end
                parent = attribute_or_na(record, "Parent")
                parent == "NA" && continue
                if code == cds_code
                    protein = attribute_or_na(record, "protein_id")
                    protein == "NA" || (cds_parents[protein] = clean(parent))
                elseif id != "NA"
                    parents[clean(id)] = clean(parent)
                end
            end
        finally
            close(reader)
        end
    end
    return (; genes, parents, cds_parents)
end

"""
Given a feature ID, follow `parents` up to a gene. Returns that gene's ID, or
`nothing` when the chain ends (or cycles) before reaching one.
"""
function _climb_to_gene(id::String, parents::Dict{String,String}, genes::Set{String})
    for _ = 0:length(parents)
        id in genes && return id
        next_id = get(parents, id, nothing)
        isnothing(next_id) && return nothing
        id = next_id
    end
    return nothing
end

"""
    protein_gene_map(gff_path; sanitize_ids = true)

Given a GFF3 file (optionally gzipped), follow each CDS's `protein_id` up its
`Parent` chain (CDS → mRNA → gene, or CDS → gene) to a `gene` record. Returns a
`protein_id => gene ID` `Dict`, IDs sanitized as [`add_features!`](@ref) would, so
they match a genome loaded with the same `sanitize_ids`. Exact where
`"gene-" * symbol` is not (`gene-TRNAK-UUU-2`); a chain reaching no gene is left out.
"""
function protein_gene_map(gff_path::String; sanitize_ids::Bool = true)
    links = _protein_links(gff_path, sanitize_ids)
    proteins = Dict{String,String}()
    for (protein, parent) in links.cds_parents
        gene = _climb_to_gene(parent, links.parents, links.genes)
        isnothing(gene) || (proteins[protein] = gene)
    end
    return proteins
end

#= Base.show overloads =#

Base.show(io::IO, f::FeatureRecord) = print(
    io,
    "FeatureRecord(\"$(f.id)\", $(f.feature_type), $(f.chromosome):$(f.start_pos)-$(f.end_pos), metadata=$(f.metadata))",
)

function Base.show(io::IO, s::Scaffold)
    n_features = length(s.features)
    print(io, "Scaffold(\"$(s.name)\", $(n_features) feature$(n_features == 1 ? "" : "s"))")
end

Base.show(io::IO, r::ParseResult) = print(
    io,
    "ParseResult($(r.scaffold_id):$(r.start_pos)-$(r.end_pos), id=\"$(r.id)\", biotype=$(r.biotype))",
)

function Base.show(io::IO, g::Genome)
    n_scaffolds = length(g.scaffolds)
    n_features = sum(length(sc.features) for sc in values(g.scaffolds); init = 0)
    print(
        io,
        "Genome($(n_scaffolds) scaffold$(n_scaffolds == 1 ? "" : "s"), $(n_features) feature$(n_features == 1 ? "" : "s"))",
    )
end

function Base.show(io::IO, sp::Species)
    taxon = isempty(sp.taxon_id) ? "" : ", taxon=$(sp.taxon_id)"
    print(io, "Species(\"$(sp.name)\"$(taxon), $(sp.genome))")
end

export FeatureRecord,
    Genome,
    IntervalTreeM64,
    ParseResult,
    Scaffold,
    Species,
    add_features!,
    convert_so_term,
    feature_id,
    feature_record,
    get_feature,
    get_metadata,
    get_metadata_id,
    get_so_terms,
    get_strand,
    pack_metadata,
    parse_index,
    parse_record,
    parse_so_term,
    parse_strand,
    protein_gene_map,
    sanitize_id,
    strand_code

end
