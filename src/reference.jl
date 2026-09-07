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
The top level of an in-memory genome representation.

- `scaffolds`: the sequences, keyed by name.
- `vocab`/`vocab_lookup`: string intern pool for feature metadata.
- `meta_offsets`/`meta_blob`: per-feature offsets into a byte blob of interned
  metadata tokens.

## Metadata handling
Each feature's metadata is split between the 64-bit code stored on its interval
and the byte blob. Bits 0-31 of the code are a 1-based index into
`meta_offsets`, which gives the feature's start offset in `meta_blob`; the
remaining bits hold the strand and SO term (see [`pack_metadata`](@ref)). A
feature's blob entry is a run of `UInt32` vocabulary tokens — ID, source and
biotype — resolved through `vocab`.
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
    while true
        matched = match(r"^(gene|transcript|mrna|rna|cds|protein):(.+)$"i, cleaned)
        matched === nothing && return cleaned
        cleaned = String(matched.captures[2])
    end
end

"""Single value of GFF3 attribute `key`, or `"NA"` when absent or multi-valued."""
function attribute_or_na(attributes::Dict, key::String)
    found = get(attributes, key, nothing)
    return (found === nothing || length(found) != 1) ? "NA" : only(found)
end

"""
Reduce a GFF3 record to a [`ParseResult`](@ref) carrying metadata index
`meta_index`, or `nothing` when its feature type is not a known SO term.
Feature IDs are stripped of type prefixes unless `sanitize_ids` is `false`.
"""
function parse_record(
    record::GFF3.Record,
    meta_index::UInt32;
    sanitize_ids::Bool = true,
)::Union{Nothing,ParseResult}
    so_term = convert_so_term(GFF3.featuretype(record))
    isnothing(so_term) && return nothing

    attributes = Dict(GFF3.attributes(record))
    feature_id = attribute_or_na(attributes, "ID")
    sanitize_ids && (feature_id = sanitize_id(feature_id))

    return ParseResult(
        GFF3.seqname(record),
        UInt32(GFF3.seqstart(record)),
        UInt32(GFF3.seqend(record)),
        pack_metadata(meta_index, strand_code(GFF3.strand(record)), so_term),
        feature_id,
        GFF3.hassource(record) ? GFF3.source(record) : "NA",
        attribute_or_na(attributes, "gene_biotype"),
    )
end

#= Metadata store =#

# Read the UInt32 vocabulary token stored at 1-based byte `offset` of `blob`.
@inline _read_token(blob::Vector{UInt8}, offset::Integer) =
    GC.@preserve blob unsafe_load(Ptr{UInt32}(pointer(blob, offset)))

# Append `token`'s 4 bytes to `blob` (native endian, matching `_read_token`).
@inline function _append_token!(blob::Vector{UInt8}, token::UInt32)
    offset = length(blob) + 1
    resize!(blob, offset + 3)
    GC.@preserve blob unsafe_store!(Ptr{UInt32}(pointer(blob, offset)), token)
    return blob
end

# Byte range of `meta_index`'s metadata entry, or `nothing` when the index is
# out of range or the entry is empty.
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

# Lazily iterate `(scaffold_name, interval)` over every feature in `genome`;
# the shared walk behind the `getindex` methods below.
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
        if get_metadata_id(genome, parse_index(interval.value)) == id
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
        feature_id = get_metadata_id(genome, parse_index(interval.value))
        if !isnothing(feature_id) && feature_id in wanted
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

# Interns `s` into the genome's vocab, returning its 1-based UInt32 token.
function _intern_string!(genome::Genome, s::String)
    get!(genome.vocab_lookup, s) do
        token = UInt32(length(genome.vocab) + 1)
        push!(genome.vocab, s)
        token
    end
end

# Runs on a dedicated task. Drains batches of ParseResults from `ch` and commits
# them into `genome` (intervals + metadata blob).
function _build_genome!(ch::Channel{Vector{ParseResult}}, genome::Genome)
    for batch in ch
        for result in batch
            scaffold = get!(genome.scaffolds, result.scaffold_id) do
                Scaffold(result.scaffold_id, IntervalTreeM64())
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
    end

    # Sentinel offset so the last feature's end byte can be computed.
    push!(genome.meta_offsets, UInt32(length(genome.meta_blob) + 1))
end

"""
Parse the GFF3 file at `gff_path` (optionally gzipped) and add every feature
whose type is a known SO term to `genome`.

The main task batches parsed records onto a channel while a spawned task commits
them, so parsing and tree construction overlap. Feature IDs are stripped of type
prefixes unless `sanitize_ids` is `false` (see [`sanitize_id`](@ref)).
"""
function add_features!(gff_path::String, genome::Genome; sanitize_ids::Bool = true)
    # Unbounded channel so the parser never blocks waiting for the builder
    ch = Channel{Vector{ParseResult}}(Inf)
    builder_task = Threads.@spawn _build_genome!(ch, genome)

    open(gff_path) do fh
        reader =
            endswith(gff_path, ".gz") ? GFF3.Reader(GzipDecompressorStream(fh)) :
            GFF3.Reader(fh)

        record = GFF3.Record()
        # meta_index is 1-based: it becomes the array index used by get_metadata
        meta_index = UInt32(1)
        buffer = sizehint!(ParseResult[], RECORD_BUFFER)

        try
            while !eof(reader)
                # No need to `empty!` record; GFF3.jl does that inside `read!`
                read!(reader, record)
                BioGenerics.isfilled(record) || continue

                result = parse_record(record, meta_index; sanitize_ids = sanitize_ids)
                isnothing(result) && continue

                push!(buffer, result)
                meta_index += UInt32(1)
                if length(buffer) == RECORD_BUFFER
                    put!(ch, buffer)
                    buffer = sizehint!(ParseResult[], RECORD_BUFFER)
                end
            end
            # Flush any remaining records that didn't fill a complete batch
            isempty(buffer) || put!(ch, buffer)
        finally
            close(reader)
            close(ch)
        end
    end

    wait(builder_task)
end

function add_features!(gff_path::String, species::Species; sanitize_ids::Bool = true)
    add_features!(gff_path, species.genome; sanitize_ids = sanitize_ids)
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
    sanitize_id,
    strand_code

end
