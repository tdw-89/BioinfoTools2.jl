module Reference

using BioGenerics
using CodecZlib
using GFF3
using IntervalTrees

using ..SOTerms

const RECORD_BUFFER = 1_000
const IntervalMeta64 = IntervalTree{UInt32, IntervalValue{UInt32, UInt64}}

##############################= Scaffolds =##############################

"""
- `name`: The name of the scaffold
- `features`: An interval tree containing the start and end of each feature, along with a 64 bit metadata code (see [below](#metadata-handling)).
"""
struct Scaffold
    # Scaffold metadata
    name::String

    # Intervals
    features::IntervalMeta64
end

function convert_strand(strand::GFF3.GenomicFeatures.Strand)
    if strand == GFF3.GenomicFeatures.STRAND_POS
        UInt8(1)
    elseif strand == GFF3.GenomicFeatures.STRAND_NEG
        UInt8(2)
    elseif strand == GFF3.GenomicFeatures.STRAND_BOTH
        UInt8(3)
    else
        UInt8(0)
    end
end

"""Get the 16-bit code for a given SO term label"""
function convert_so_term(label::String)
    so_term_result = SO_TERMS[label]
    if isnothing(so_term_result)
        return nothing
    end 
    return so_term_result[1]
end

"""
|-<NULL>-|-----SO-Term-----|-Strand-|---------------Index---------------|
|00000000|00000000|00000000|00000000|00000000|00000000|00000000|00000000|
"""
function pack_metadata(index::UInt32, strand::UInt8, so_term::UInt16)
    code = UInt64(0)
    code = code | UInt64(index)
    code = code | (UInt64(strand) << 32)
    code = code | (UInt64(so_term) << (32 + 8))
    return code
end

function parse_index(code::UInt64)
    return UInt32(code & 0x00000000FFFFFFFF)
end

function parse_strand(code::UInt64)
    bit_code = UInt8((code >> 32) & 0xFF)
    if bit_code == 1
        return GFF3.GenomicFeatures.STRAND_POS
    elseif bit_code == 2
        return GFF3.GenomicFeatures.STRAND_NEG
    elseif bit_code == 3
        return GFF3.GenomicFeatures.STRAND_BOTH
    else
        return GFF3.GenomicFeatures.STRAND_NA
    end
end

function parse_so_term(code::UInt64)
    return UInt16(code >> (32 + 8)) & 0xFFFF
end

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
Strip common feature-type prefixes (for example `gene:` or `transcript:`)
from IDs when present.
"""
function sanitize_id(id::AbstractString)
    cleaned = String(id)
    # Some annotations can stack prefixes (e.g. "gene:transcript:...").
    while true
        m = match(r"^(gene|transcript|mrna|rna|cds|protein):(.+)$"i, cleaned)
        m === nothing && break
        cleaned = m.captures[2]
    end
    return cleaned
end

function parse_record(record::GFF3.Record, meta_index::UInt32; sanitize_ids::Bool = true)
    # interned
    scaffold = GFF3.seqname(record)
    feature_attr = GFF3.attributes(record) |> Dict
    feature_id = haskey(feature_attr, "ID") && length(feature_attr["ID"]) == 1 ? only(feature_attr["ID"]) : "NA"
    if sanitize_ids
        feature_id = sanitize_id(feature_id)
    end
    feature_source = GFF3.hassource(record) ? GFF3.source(record) : "NA"
    gene_biotype = haskey(feature_attr, "gene_biotype") && length(feature_attr["gene_biotype"]) == 1 ? only(feature_attr["gene_biotype"]) : "NA"
    
    # bits
    feature_type = record |> GFF3.featuretype
    so_term = feature_type |> convert_so_term
    if isnothing(so_term)
        return nothing
    end
    strand = record |> GFF3.strand |> convert_strand
    code = pack_metadata(meta_index, strand, so_term)
    
    # Interval
    start_pos = record |> GFF3.seqstart |> UInt32
    end_pos = record |> GFF3.seqend |> UInt32

    return ParseResult(
        scaffold,
        start_pos,
        end_pos,
        code,
        feature_id,
        feature_source,
        gene_biotype
    )
end

"""
This type is the top level of organization 
for an in-memory repressentation of a genome.

At the time of writing it has the following components:
- `vocab`/`vocab_lookup`: The vocabulary and lookup for string interning extended metadata
- `meta_offsets`/`meta_blobs`: The stored offsets and actual metadata byte blobs

## Metadata handling
The metadata for each feature is stored using the 64-bit 'code' and the byte ([`UInt8`](@ref)) blob.
The first 32 bits contains parseable metadata (currently strand and SO term, 8 bits still unused)
while the second 32 bits are an index into the `meta_offsets` vector which gives the offset into the
byte blob for that feature.
"""
mutable struct Genome
    scaffolds::Dict{String, Scaffold}

    # String intern pool
    vocab::Vector{String}
    vocab_lookup::Dict{String, UInt32}

    # Metadata store
    meta_offsets::Vector{UInt32}
    meta_blob::Vector{UInt8}
end

function get_metadata(genome::Genome, meta_index::UInt32)
    if length(genome.meta_offsets) <= 1 || (meta_index + 1) > length(genome.meta_offsets)
        return String[]
    end

    start_byte = genome.meta_offsets[meta_index]
    end_byte = genome.meta_offsets[meta_index + 1] - 1
    
    if start_byte > end_byte
        return String[]
    end
    
    raw_bytes = genome.meta_blob[start_byte:end_byte]
    tokens = reinterpret(UInt32, raw_bytes)
    
    return [genome.vocab[t] for t in tokens]
end

"""
Return only the interned ID (the first metadata token) for `meta_index`, or
`nothing` when the feature carries no metadata.

The ID is always stored as the first `UInt32` token of a feature's metadata
entry, so this reads just those 4 bytes directly from the blob (native endian,
matching how `build_genome!` wrote it) instead of allocating the full metadata
vector. Prefer this over `get_metadata` when only the ID is needed.
"""
function get_metadata_id(genome::Genome, meta_index::UInt32)
    if meta_index < 0x1 || length(genome.meta_offsets) <= 1 || (meta_index + 1) > length(genome.meta_offsets)
        return nothing
    end

    start_byte = genome.meta_offsets[meta_index]
    end_byte = genome.meta_offsets[meta_index + 1] - 1

    if start_byte > end_byte
        return nothing
    end

    token = GC.@preserve genome unsafe_load(Ptr{UInt32}(pointer(genome.meta_blob, start_byte)))
    return genome.vocab[token]
end

function get_metadata(genome::Genome, features::IntervalMeta64)
    metadata_list = Vector{Vector{String}}(undef, length(features))
    for (i, interval) in enumerate(features)
        meta_idx = parse_index(interval.value)
        metadata_list[i] = get_metadata(genome, meta_idx)
    end
    return metadata_list
end

get_metadata(genome::Genome, scaffold::Scaffold) = get_metadata(genome, scaffold.features)

function get_metadata(genome::Genome)
    scaffolds = Dict{String, Vector{Vector{String}}}()
    for (scaffold_name, scaffold) in genome.scaffolds
        scaffolds[scaffold_name] = get_metadata(genome, scaffold)
    end
    return scaffolds
end

function get_feature(scaffold::Scaffold, feature::Symbol)::Union{Nothing, IntervalMeta64}
    result = SO_TERMS[feature]
    isnothing(result) && return nothing
    feature_bit_mask, _ = result

    tree = IntervalMeta64()
    for interval in scaffold.features
        if parse_so_term(interval.value) == feature_bit_mask
            push!(tree, interval)
        end
    end
    return tree
end

get_feature(scaffold::Scaffold, feature::AbstractString)::Union{Nothing, IntervalMeta64} = get_feature(scaffold, Symbol(feature))

function get_feature(genome::Genome, feature::Symbol)
    result = SO_TERMS[feature]
    isnothing(result) && return Dict{String, IntervalMeta64}()

    scaffolds = Dict{String, IntervalMeta64}()
    for (scaffold_name, scaffold) in genome.scaffolds
        scaffolds[scaffold_name] = get_feature(scaffold, feature)
    end
    return scaffolds
end

get_feature(genome::Genome, feature::AbstractString) = get_feature(genome, Symbol(feature))

"""
Return the unique SO term labels used by features stored in `genome`.
"""
function get_so_terms(genome::Genome)
    terms = Set{Symbol}()
    for scaffold in values(genome.scaffolds)
        for interval in scaffold.features
            short_id = parse_so_term(interval.value)
            term = SO_TERMS[short_id]
            if !isnothing(term)
                push!(terms, term[2])
            end
        end
    end
    return sort!(collect(terms))
end

"""
A single annotated feature looked up from a `Genome`:
- `id`: the ID the feature was searched by (first metadata field)
- `feature_type`: the SO term label (e.g. `:gene`)
- `chromosome`: the scaffold name the feature lives on
- `start_pos` / `end_pos`: 1-based closed interval bounds
- `metadata`: the parsed metadata stored for the feature (ID, source, biotype)
"""
struct FeatureRecord
    id::String
    feature_type::Symbol
    chromosome::String
    start_pos::UInt32
    end_pos::UInt32
    metadata::Vector{String}
end

# Build a FeatureRecord from a scaffold name and one of its interval entries.
function _feature_record(genome::Genome, chromosome::String, interval::IntervalValue{UInt32, UInt64})
    code = interval.value
    so_result = SO_TERMS[parse_so_term(code)]
    feature_type = isnothing(so_result) ? Symbol("") : so_result[2]
    metadata = get_metadata(genome, parse_index(code))
    id = isempty(metadata) ? "" : metadata[1]
    return FeatureRecord(id, feature_type, chromosome, interval.first, interval.last, metadata)
end

"""
Look up a single feature by its ID (the first metadata field). Every feature's
ID is resolved on demand and compared to `id`; the first match is returned as a
[`FeatureRecord`](@ref), or `nothing` when no feature matches.

This scans the genome's features (O(features) per lookup) rather than keeping an
ID index.
"""
function Base.getindex(genome::Genome, id::AbstractString)
    for (name, scaffold) in genome.scaffolds
        for interval in scaffold.features
            if get_metadata_id(genome, parse_index(interval.value)) == id
                return _feature_record(genome, name, interval)
            end
        end
    end
    return nothing
end

"""
Look up every feature whose ID (first metadata field) is contained in `ids`,
returning a `Vector{FeatureRecord}`. IDs with no matching feature are skipped.
"""
function Base.getindex(genome::Genome, ids::AbstractVector{<:AbstractString})
    wanted = Set{String}(ids)
    records = FeatureRecord[]
    for (name, scaffold) in genome.scaffolds
        for interval in scaffold.features
            fid = get_metadata_id(genome, parse_index(interval.value))
            if !isnothing(fid) && fid in wanted
                push!(records, _feature_record(genome, name, interval))
            end
        end
    end
    return records
end

"""
Look up the feature carrying the 32-bit metadata index `meta_index`, returning
it as a [`FeatureRecord`](@ref), or `nothing` when no feature uses that index.

Scans the genome's intervals (O(features) per lookup), comparing each interval's
parsed metadata index (see [`parse_index`](@ref)) against `meta_index`. This is a
convenience for turning the metadata index stored on a `TabularData` sample back
into a full feature record.
"""
function Base.getindex(genome::Genome, meta_index::UInt32)
    for (name, scaffold) in genome.scaffolds
        for interval in scaffold.features
            if parse_index(interval.value) == meta_index
                return _feature_record(genome, name, interval)
            end
        end
    end
    return nothing
end

function Base.show(io::IO, f::FeatureRecord)
    print(io, "FeatureRecord(\"$(f.id)\", $(f.feature_type), $(f.chromosome):$(f.start_pos)-$(f.end_pos), metadata=$(f.metadata))")
end

mutable struct Species
    name::String
    taxon_id::String
    genome::Genome
end

function Species(name::String; taxon_id::String = "")
    genome = Genome(
        Dict{String, Scaffold}(),
        String[],
        Dict{String, UInt32}(),
        UInt32[],
        UInt8[]
    )
    return Species(name, taxon_id, genome)
end

# Interns `s` into the genome's vocab, returning its 1-based UInt32 token.
function intern_string!(genome::Genome, s::String)
    get!(genome.vocab_lookup, s) do
        token = UInt32(length(genome.vocab) + 1)
        push!(genome.vocab, s)
        token
    end
end

# Runs on a dedicated CPU thread. Drains batches of ParseResults from `ch`
# and commits them into `genome` (intervals + metadata blob).
function build_genome!(ch::Channel{Vector{ParseResult}}, genome::Genome)
    for batch in ch
        for result in batch
            # Ensure the scaffold exists, creating it lazily if not
            scaffold = get!(genome.scaffolds, result.scaffold_id) do
                Scaffold(result.scaffold_id, IntervalMeta64())
            end

            # 1. Add the feature interval (start, end, 64-bit code) to the scaffold tree
            push!(scaffold.features, IntervalValue(result.start_pos, result.end_pos, result.code))

            # 2. Record the start offset for this features's metadata entry, then encode it.
            #    meta_offsets grows one entry per feature; a final sentinel is appended at the
            #    end so get_metadata can compute end_byte = meta_offsets[i+1] - 1.
            push!(genome.meta_offsets, UInt32(length(genome.meta_blob) + 1))
            for s in (result.id, result.source, result.biotype)
                token = intern_string!(genome, s)
                # Write the UInt32 token as 4 bytes (native endian, matches reinterpret in get_metadata)
                append!(genome.meta_blob, reinterpret(UInt8, [token]))
            end
        end
    end

    # Sentinel offset so the last feature's end byte can be computed by get_metadata
    push!(genome.meta_offsets, UInt32(length(genome.meta_blob) + 1))
end

function add_features!(gff_path::String, genome::Genome; sanitize_ids::Bool = true)
    # Unbounded channel so the parser never blocks waiting for the builder
    ch = Channel{Vector{ParseResult}}(Inf)
    builder_task = Threads.@spawn build_genome!(ch, genome)

    open(gff_path) do fh
        rdr = endswith(gff_path, ".gz") ?
            GFF3.Reader(GzipDecompressorStream(fh)) :
            GFF3.Reader(fh)

        record = GFF3.Record()
        # meta_index is 1-based: it becomes the array index used by get_metadata
        meta_index = UInt32(1)
        buffer = Vector{ParseResult}()
        sizehint!(buffer, RECORD_BUFFER)

        try
            while !eof(rdr)
                # Don't need to `empty!` record because GFF3.jl does this
                # as a first step in `read!`
                read!(rdr, record)
                if BioGenerics.isfilled(record)
                    result = parse_record(record, meta_index; sanitize_ids = sanitize_ids)
                    if !isnothing(result) 
                        push!(buffer, result)
                        meta_index += UInt32(1)
                        if length(buffer) == RECORD_BUFFER
                            put!(ch, buffer)
                            buffer = Vector{ParseResult}()
                            sizehint!(buffer, RECORD_BUFFER)
                        end
                    end
                end
                # at_eof && break
            end
            # Flush any remaining records that didn't fill a complete batch
            isempty(buffer) || put!(ch, buffer)
        finally
            close(rdr)
            close(ch)
        end
    end

    wait(builder_task)
end

function add_features!(gff_path::String, species::Species; sanitize_ids::Bool = true)
    add_features!(gff_path, species.genome; sanitize_ids = sanitize_ids)
end

Base.show(io::IO, s::Scaffold) = print(io, "Scaffold(\"$(s.name)\", $(length(s.features)) feature$(length(s.features) == 1 ? "" : "s"))")

function Base.show(io::IO, r::ParseResult)
    print(io, "ParseResult($(r.scaffold_id):$(r.start_pos)-$(r.end_pos), id=\"$(r.id)\", biotype=$(r.biotype))")
end

function Base.show(io::IO, g::Genome)
    nscaff = length(g.scaffolds)
    nfeatures = sum(length(sc.features) for sc in values(g.scaffolds); init=0)
    print(io, "Genome($(nscaff) scaffold$(nscaff == 1 ? "" : "s"), $(nfeatures) feature$(nfeatures == 1 ? "" : "s"))")
end

function Base.show(io::IO, sp::Species)
    taxon = isempty(sp.taxon_id) ? "" : ", taxon=$(sp.taxon_id)"
    print(io, "Species(\"$(sp.name)\"$(taxon), $(sp.genome))")
end

export
    Species,
    IntervalMeta64,
    Genome,
    Scaffold,
    FeatureRecord,
    add_features!,
    get_metadata,
    get_feature,
    get_so_terms
    
end