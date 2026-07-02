using BioinfoTools2.Studies
using Test

const ST_DATA_DIR = joinpath(@__DIR__, "data")
const MICRO_BED   = joinpath(ST_DATA_DIR, "micro.narrowPeak")
const FULL_BED    = joinpath(ST_DATA_DIR, "full.narrowPeak")

@testset "Studies" begin

    # -------------------------------------------------------------------------
    @testset "pack_bed_code / parse_bed_strand roundtrip" begin
        @test Studies.parse_bed_strand(Studies.pack_bed_code(UInt8(0))) == UInt8(0)
        @test Studies.parse_bed_strand(Studies.pack_bed_code(UInt8(1))) == UInt8(1)
        @test Studies.parse_bed_strand(Studies.pack_bed_code(UInt8(2))) == UInt8(2)

        # Only bits 33-40 should be set; all other bits must remain zero
        @test (Studies.pack_bed_code(UInt8(1)) & ~(UInt64(0xFF) << 32)) == UInt64(0)
        @test (Studies.pack_bed_code(UInt8(2)) & ~(UInt64(0xFF) << 32)) == UInt64(0)
    end

    # -------------------------------------------------------------------------
    @testset "load_bed - micro.narrowPeak (8 records, 1 scaffold)" begin
        bd = Studies.load_bed(MICRO_BED)

        @test bd isa Studies.BedData
        @test length(bd.scaffolds) == 1
        @test haskey(bd.scaffolds, "DDB0215018")
        @test length(bd.scaffolds["DDB0215018"]) == 8
    end

    # -------------------------------------------------------------------------
    @testset "load_bed - strand encoding in micro.narrowPeak" begin
        # All 8 records in micro.narrowPeak have strand '.' → UInt8(0)
        bd = Studies.load_bed(MICRO_BED)

        for (_, tree) in bd.scaffolds
            for iv in tree
                @test Studies.parse_bed_strand(iv.value) == UInt8(0)
            end
        end
    end

    # -------------------------------------------------------------------------
    @testset "load_bed - coordinate conversion in micro.narrowPeak" begin
        # BED is 0-based half-open; loader converts to 1-based closed.
        # First record: chromStart=4, chromEnd=1609 → start=5, end=1609
        bd = Studies.load_bed(MICRO_BED)
        tree = bd.scaffolds["DDB0215018"]

        starts = sort([iv.first for iv in tree])
        ends   = sort([iv.last  for iv in tree])

        @test minimum(starts) == UInt32(5)    # 0-based 4  → 1-based 5
        @test ends[1]         == UInt32(1609) # chromEnd unchanged
    end

    # -------------------------------------------------------------------------
    @testset "load_bed - full.narrowPeak (4448 records total)" begin
        bd = Studies.load_bed(FULL_BED)

        @test bd isa Studies.BedData
        @test !isempty(bd.scaffolds)

        total = sum(length(tree) for tree in values(bd.scaffolds))
        @test total == 4448
    end

end
