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
- `genes`: An interval tree containing the TSS and TES of each gene, along with a 64 bit metadata code (see [below](#metadata-handling)).
"""
struct Scaffold
    # Scaffold metadata
    name::String

    # Intervals
    genes::IntervalMeta64
end


function convert_strand(strand::GFF3.GenomicFeatures.Strand)
    if strand == GFF3.GenomicFeatures.STRAND_NA
        UInt8(0)
    elseif strand == GFF3.GenomicFeatures.STRAND_POS
        UInt8(1)
    elseif strand == GFF3.GenomicFeatures.STRAND_NEG
        UInt8(2)
    elseif strand == GFF3.GenomicFeatures.STRAND_NEG
        UInt8(3)
    end
end

"""Get the 16-bit code for a given SO term label"""
function convert_so_term(label::String)
    bit_code, _ = SO_TERMS[label]
    return bit_code
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
    return UInt8((code >> 32) & 0xFF)
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

"""WIP: Currently only 'gene' SO type supported"""
function parse_record(record::GFF3.Record, meta_index::UInt32)
    if GFF3.featuretype(record) == "gene"
        # interned
        scaffold = GFF3.seqname(record)
        gene_attr = GFF3.attributes(record) |> Dict
        gene_id = haskey(gene_attr, "ID") && length(gene_attr["ID"]) == 1 ? only(gene_attr["ID"]) : "NA"
        gene_source = GFF3.source(record)
        gene_biotype = haskey(gene_attr, "gene_biotype") && length(gene_attr["gene_biotype"]) == 1 ? only(gene_attr["gene_biotype"]) : "NA"
        
        # bits
        so_term = record |> GFF3.featuretype |> convert_so_term
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
            gene_id,
            gene_source,
            gene_biotype
        )
    end
end

"""
This type is the top level of organization 
for an in-memory repressentation of a genome.

At the time of writing it has the following components:
- `vocab`/`vocab_lookup`: The vocabulary and lookup for string interning extended metadata
- `meta_offsets`/`meta_blobs`: The stored offsets and actual metadata byte blobs

## Metadata handling
The metadata for each gene is stored using the 64-bit 'code' and the byte ([`UInt8`](@ref)) blob.
The first 32 bits contains parseable metadata (currently strand only, but soon other items like biotype)
while the second 32 bits are an index into the `meta_offsets` vector which gives the offset into the
byte blob for that gene.
"""
mutable struct Genome
    scaffolds::Dict{String, Scaffold}

    # String intern pool
    vocab::Vector{String}
    vocab_lookup::Dict{String, UInt16}

    # Metadata store
    meta_offsets::Vector{UInt32}
    meta_blob::Vector{UInt8}
end


function get_metadata(genome::Genome, meta_index::UInt32)
    # 1. Look up the byte boundaries
    start_byte = genome.meta_offsets[meta_index]
    end_byte = genome.meta_offsets[meta_index + 1] - 1
    
    # If there are no tags for this gene, return an empty array
    if start_byte > end_byte
        return String[]
    end
    
    # 2. Slice the raw bytes (this creates a fast, contiguous copy)
    raw_bytes = genome.meta_blob[start_byte:end_byte]
    
    # 3. Reinterpret the bytes as UInt16 tokens
    tokens = reinterpret(UInt16, raw_bytes)
    
    # 4. Map the tokens back to human-readable strings
    return [genome.vocab[t] for t in tokens]
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
        Dict{String, UInt16}(),
        UInt32[],
        UInt8[]
    )
    return Species(name, taxon_id, genome)
end

# Interns `s` into the genome's vocab, returning its 1-based UInt16 token.
function intern_string!(genome::Genome, s::String)
    get!(genome.vocab_lookup, s) do
        token = UInt16(length(genome.vocab) + 1)
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

            # 1. Add the gene interval (start, end, 64-bit code) to the scaffold tree
            push!(scaffold.genes, IntervalValue(result.start_pos, result.end_pos, result.code))

            # 2. Record the start offset for this gene's metadata entry, then encode it.
            #    meta_offsets grows one entry per gene; a final sentinel is appended at the
            #    end so get_metadata can compute end_byte = meta_offsets[i+1] - 1.
            push!(genome.meta_offsets, UInt32(length(genome.meta_blob) + 1))
            for s in (result.id, result.source, result.biotype)
                token = intern_string!(genome, s)
                # Write the UInt16 token as 2 bytes (native endian, matches reinterpret in get_metadata)
                append!(genome.meta_blob, reinterpret(UInt8, [token]))
            end
        end
    end

    # Sentinel offset so the last gene's end byte can be computed by get_metadata
    push!(genome.meta_offsets, UInt32(length(genome.meta_blob) + 1))
end

function add_features!(gff_path::String, genome::Genome)
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
                    result = parse_record(record, meta_index)
                    if result !== nothing
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

function add_features!(gff_path::String, species::Species)
    add_features!(gff_path, species.genome)
end

export
    Species,
    Genome,
    Scaffold,
    add_features!,
    get_metadata
    
end