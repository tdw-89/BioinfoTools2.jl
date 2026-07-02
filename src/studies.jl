module Studies

using BED
using BioGenerics
using CodecZlib
using Dates
using GFF3
using IntervalTrees
using ..Reference

mutable struct BedData
    scaffolds::Dict{String, Reference.IntervalMeta64}

end

"""
Bit layout for the 64-bit metadata code attached to each BED interval.
Bits 33-40 encode strand using the same scheme as gene intervals in `Scaffold`.

|  64-57  |  56-41  | 40-33  |          32-1          |
|---------|---------|--------|------------------------|
| <NULL>  | <NULL>  | Strand |        <NULL>          |

- Bits  1-32  : reserved
- Bits 33-40  : strand (0 = unknown/unstranded, 1 = +, 2 = -)
- Bits 41-64  : reserved
"""
function pack_bed_code(strand::UInt8)::UInt64
    return UInt64(strand) << 32
end

function parse_bed_strand(code::UInt64)::UInt8
    return UInt8((code >> 32) & 0xFF)
end

"""
Convert a BED record's strand field to the same UInt8 encoding used by genes.
"""
function bed_record_strand(record::BED.Record)
    BED.hasstrand(record) || return UInt8(0)
    s = BED.strand(record)
    if s == GFF3.GenomicFeatures.STRAND_POS
        UInt8(1)
    elseif s == GFF3.GenomicFeatures.STRAND_NEG
        UInt8(2)
    else
        UInt8(0)
    end
end

function load_bed(file_path::String)
    scaffolds = Dict{String, Reference.IntervalMeta64}()

    open(file_path) do fh
        rdr = endswith(file_path, ".gz") ?
            BED.Reader(GzipDecompressorStream(fh)) :
            BED.Reader(fh)

        record = BED.Record()
        try
            while true
                empty!(record)
                at_eof = false
                try
                    read!(rdr, record)
                catch e
                    e isa EOFError || rethrow()
                    at_eof = true
                end
                # Process if filled: handles both normal reads and the last record
                # in files without a trailing newline (EOFError thrown after fill).
                if BED.isfilled(record)
                    # BED uses 0-based half-open [start, end) coordinates;
                    # convert to 1-based closed [start, end] to match gene intervals.
                    # from GFF3 files.
                    chrom     = BED.chrom(record)
                    start_pos = UInt32(BED.chromstart(record) + 1)
                    end_pos   = UInt32(BED.chromend(record))

                    strand = bed_record_strand(record)
                    code   = pack_bed_code(strand)

                    tree = get!(Reference.IntervalMeta64, scaffolds, chrom)
                    push!(tree, IntervalValue(start_pos, end_pos, code))
                end
                at_eof && break
            end
        finally
            close(rdr)
        end
    end

    return BedData(scaffolds)
end

mutable struct TabularData{T}
    table::Matrix{T}
end

const Data = Union{BedData, TabularData}

mutable struct BioSample
    sample_id::String
    tissue_type::String
    species::Species
end

mutable struct AssayMethod
    name::String
    description::String
end

mutable struct Measurement
    file_path::String
    format::String
    data::Data
end

mutable struct Assay
    id::String
    type::String
    description::String
    measurement::Measurement
    biosample::BioSample
    method::AssayMethod
end

mutable struct Study
    id::String
    title::String
    date::Date
    assays::Vector{Assay}
end

end