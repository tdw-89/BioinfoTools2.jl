using BioinfoTools2.Studies
using BioinfoTools2.Reference
using IntervalTrees
using Test

const ST_DATA_DIR      = joinpath(@__DIR__, "data")
const MICRO_BED        = joinpath(ST_DATA_DIR, "micro.bed")
const MICRO_NARROWPEAK = joinpath(ST_DATA_DIR, "micro.narrowPeak")
const FULL_NARROWPEAK  = joinpath(ST_DATA_DIR, "full.narrowPeak")
const ST_GFF_SINGLE    = joinpath(ST_DATA_DIR, "NC_003280.10.gff.gz")

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
        bd = Studies.load_bed(MICRO_NARROWPEAK)

        @test bd isa Studies.BedData
        @test length(bd.scaffolds) == 1
        @test haskey(bd.scaffolds, "DDB0215018")
        @test length(bd.scaffolds["DDB0215018"]) == 8
    end

    # -------------------------------------------------------------------------
    @testset "load_bed - strand encoding in micro.narrowPeak" begin
        # All 8 records in micro.narrowPeak have strand '.' → UInt8(0)
        bd = Studies.load_bed(MICRO_NARROWPEAK)

        for (_, tree) in bd.scaffolds
            for iv in tree
                @test Studies.parse_bed_strand(iv.value) == UInt8(0)
            end
        end
    end

    # -------------------------------------------------------------------------
    @testset "load_bed - coordinate conversion in micro.narrowPeak" begin
        # BED is 0-based half-open; BED.jl should convert to 1-based closed.
        # First record: chromStart=4, chromEnd=1609 → start=5, end=1609
        bd = Studies.load_bed(MICRO_NARROWPEAK)
        tree = bd.scaffolds["DDB0215018"]

        starts = sort([iv.first for iv in tree])
        ends   = sort([iv.last  for iv in tree])

        @test minimum(starts) == UInt32(5)    # 0-based 4  → 1-based 5
        @test ends[1]         == UInt32(1609) # chromEnd unchanged
    end

    # -------------------------------------------------------------------------
    @testset "load_bed - BED (3 + 3)" begin
        bd = Studies.load_bed(MICRO_BED)
        tree = bd.scaffolds["DDB0215018"]

        starts = sort([iv.first for iv in tree])
        ends   = sort([iv.last  for iv in tree])

        @test minimum(starts) == UInt32(5)    # 0-based 4  → 1-based 5
        @test ends[1]         == UInt32(1609) # chromEnd unchanged
    end

    # -------------------------------------------------------------------------
    @testset "load_bed - full.narrowPeak (4448 records total)" begin
        bd = Studies.load_bed(FULL_NARROWPEAK)

        @test bd isa Studies.BedData
        @test !isempty(bd.scaffolds)

        total = sum(length(tree) for tree in values(bd.scaffolds))
        @test total == 4448
    end

    # =========================================================================
    @testset "intersect – all overloads" begin
        # Load a real genome (78 407 features on scaffold \"NC_003280.10\") once
        # for all Scaffold / Genome -based sub-tests.
        sp = Species("C. elegans")
        add_features!(ST_GFF_SINGLE, sp.genome)
        scaffold_nc = sp.genome.scaffolds["NC_003280.10"]

        # BedData whose only scaffold name matches the loaded GFF scaffold.
        bed_nc = let
            tree = IntervalMeta64()
            push!(tree, IntervalValue(UInt32(1), UInt32(100_000_000), UInt64(0)))
            Studies.BedData(Dict("NC_003280.10" => tree))
        end

        # BedData with a scaffold name that has no match in sp.genome.
        bed_no_match = Studies.BedData(Dict("NOMATCH" => IntervalMeta64()))

        # -----------------------------------------------------------------------
        @testset "intersect(IntervalMeta64, IntervalMeta64) – overlap" begin
            tree_a = IntervalMeta64()
            push!(tree_a, IntervalValue(UInt32(10),  UInt32(50),  UInt64(7)))
            push!(tree_a, IntervalValue(UInt32(200), UInt32(300), UInt64(0)))
            tree_b = IntervalMeta64()
            push!(tree_b, IntervalValue(UInt32(30), UInt32(100), UInt64(0)))

            result = Studies.intersect(tree_a, tree_b)
            ivs = collect(result)
            @test length(ivs) == 1
            @test ivs[1].first == UInt32(30)
            @test ivs[1].last  == UInt32(50)
            @test ivs[1].value == UInt64(7)   # value carried from tree_a
        end

        @testset "intersect(IntervalMeta64, IntervalMeta64) – no overlap" begin
            tree_a = IntervalMeta64()
            push!(tree_a, IntervalValue(UInt32(1), UInt32(10), UInt64(0)))
            tree_b = IntervalMeta64()
            push!(tree_b, IntervalValue(UInt32(20), UInt32(30), UInt64(0)))
            @test isempty(Studies.intersect(tree_a, tree_b))
        end

        @testset "intersect(IntervalMeta64, IntervalMeta64) – both empty" begin
            @test isempty(Studies.intersect(IntervalMeta64(), IntervalMeta64()))
        end

        # -----------------------------------------------------------------------
        @testset "intersect(Scaffold, BedData) – scaffold not in BedData" begin
            @test isnothing(Studies.intersect(scaffold_nc, bed_no_match))
        end

        @testset "intersect(Scaffold, BedData) – scaffold matches" begin
            result = Studies.intersect(scaffold_nc, bed_nc)
            @test !isnothing(result)
            @test !isempty(result)
        end

        # -----------------------------------------------------------------------
        @testset "intersect(Genome, BedData) – no matching scaffold" begin
            result = Studies.intersect(sp.genome, bed_no_match)
            @test result isa Dict
            @test isempty(result)
        end

        @testset "intersect(Genome, BedData) – matching scaffold" begin
            result = Studies.intersect(sp.genome, bed_nc)
            @test result isa Dict
            @test haskey(result, "NC_003280.10")
            @test !isempty(result["NC_003280.10"])
        end

        # -----------------------------------------------------------------------
        @testset "intersect(Scaffold, BedData, feature) – scaffold not in BedData" begin
            @test isnothing(Studies.intersect(scaffold_nc, bed_no_match, :gene))
        end

        @testset "intersect(Scaffold, BedData, feature) – valid feature, matching scaffold" begin
            result = Studies.intersect(scaffold_nc, bed_nc, :gene)
            @test !isnothing(result)
            @test !isempty(result)
        end

        @testset "intersect(Scaffold, BedData, feature) – unknown feature" begin
            @test isnothing(Studies.intersect(scaffold_nc, bed_nc, :not_a_real_so_term))
        end

        # -----------------------------------------------------------------------
        @testset "intersect(Genome, BedData, feature) – no matching scaffold" begin
            result = Studies.intersect(sp.genome, bed_no_match, :gene)
            @test result isa Dict
            @test isempty(result)
        end

        @testset "intersect(Genome, BedData, feature) – valid feature, matching scaffold" begin
            result = Studies.intersect(sp.genome, bed_nc, :gene)
            @test result isa Dict
            @test haskey(result, "NC_003280.10")
            @test !isempty(result["NC_003280.10"])
        end

        @testset "intersect(Genome, BedData, feature) – unknown feature" begin
            result = Studies.intersect(sp.genome, bed_nc, :not_a_real_so_term)
            @test result isa Dict
            @test isempty(result)
        end

    end  # intersect – all overloads
end
