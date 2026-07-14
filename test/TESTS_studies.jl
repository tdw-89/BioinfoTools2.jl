using BioinfoTools2.Studies
using BioinfoTools2.Reference
using DataFrames
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
        bd = load_bed(MICRO_NARROWPEAK)

        @test bd isa BedData
        @test length(bd.scaffolds) == 1
        @test haskey(bd.scaffolds, "DDB0215018")
        @test length(bd.scaffolds["DDB0215018"]) == 8
    end

    # -------------------------------------------------------------------------
    @testset "load_bed - strand encoding in micro.narrowPeak" begin
        # All 8 records in micro.narrowPeak have strand '.' → UInt8(0)
        bd = load_bed(MICRO_NARROWPEAK)

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
        bd = load_bed(MICRO_NARROWPEAK)
        tree = bd.scaffolds["DDB0215018"]

        starts = sort([iv.first for iv in tree])
        ends   = sort([iv.last  for iv in tree])

        @test minimum(starts) == UInt32(5)    # 0-based 4  → 1-based 5
        @test ends[1]         == UInt32(1609) # chromEnd unchanged
    end

    # -------------------------------------------------------------------------
    @testset "load_bed - BED (3 + 3)" begin
        bd = load_bed(MICRO_BED)
        tree = bd.scaffolds["DDB0215018"]

        starts = sort([iv.first for iv in tree])
        ends   = sort([iv.last  for iv in tree])

        @test minimum(starts) == UInt32(5)    # 0-based 4  → 1-based 5
        @test ends[1]         == UInt32(1609) # chromEnd unchanged
    end

    # -------------------------------------------------------------------------
    @testset "load_bed - full.narrowPeak (4448 records total)" begin
        bd = load_bed(FULL_NARROWPEAK)

        @test bd isa BedData
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
            BedData(Dict("NC_003280.10" => tree))
        end

        # BedData with a scaffold name that has no match in sp.genome.
        bed_no_match = BedData(Dict("NOMATCH" => IntervalMeta64()))

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

    # =========================================================================
    @testset "TabularData indexing & DataFrame conversion" begin
        # Build a real TabularData: load a genome, grab a few real gene IDs, then
        # assemble a DataFrame whose first column carries those IDs (plus one that
        # matches nothing) and two numeric variable columns.
        sp = Species("C. elegans")
        add_features!(ST_GFF_SINGLE, sp.genome)
        scaffold_nc = sp.genome.scaffolds["NC_003280.10"]

        gene_ids = String[]
        for iv in get_feature(scaffold_nc, :gene)
            gid = Reference.get_metadata_id(sp.genome, Reference.parse_index(iv.value))
            gid !== nothing && push!(gene_ids, gid)
            length(gene_ids) >= 3 && break
        end

        df = DataFrame(
            sample = vcat(gene_ids, ["NO_SUCH_ID"]),
            v1 = [10.0, 20.0, 30.0, 40.0],
            v2 = [11.0, 21.0, 31.0, 41.0],
        )
        tab = load_table(sp.genome, df, :gene)

        # --- Integer / cartesian indexing ------------------------------------
        @testset "getindex(Integer...) – scalar elements" begin
            @test tab[1, 1] == 10.0
            @test tab[3, 2] == 31.0
            @test tab[4, 1] == 40.0
        end

        @testset "getindex(Integer...) – row ranges and column slices" begin
            @test tab[1:2, :] == [10.0 11.0; 20.0 21.0]
            @test tab[:, 1]   == [10.0, 20.0, 30.0, 40.0]
            @test tab[4, :]   == [40.0, 41.0]
        end

        # --- Single-string lookup --------------------------------------------
        @testset "getindex(genome, String) – match returns 1-row TabularData" begin
            one = tab[sp.genome, gene_ids[2]]
            @test one isa TabularData
            @test size(one.table) == (1, 2)
            @test one.table == [20.0 21.0]
            @test one.variables == tab.variables
            @test length(one.samples) == 1
            @test Studies.sample_id(sp.genome, one.samples[1]) == gene_ids[2]
        end

        @testset "getindex(genome, String) – no match returns nothing" begin
            # "NO_SUCH_ID" is present in the table but never matched a feature,
            # so its sample slot is `nothing` and cannot be looked up.
            @test isnothing(tab[sp.genome, "NO_SUCH_ID"])
            @test isnothing(tab[sp.genome, "totally_absent_id"])
        end

        # --- Vector-of-strings lookup ----------------------------------------
        @testset "getindex(genome, Vector{String}) – subset, skipping misses" begin
            # Query order is deliberately scrambled and includes an absent ID;
            # the result must follow the original table's row order.
            sub = tab[sp.genome, [gene_ids[3], "absent_id", gene_ids[1]]]
            @test sub isa TabularData
            @test size(sub.table) == (2, 2)
            @test sub.table == [10.0 11.0; 30.0 31.0]
            @test sub.variables == tab.variables
            @test Studies.sample_id(sp.genome, sub.samples[1]) == gene_ids[1]
            @test Studies.sample_id(sp.genome, sub.samples[2]) == gene_ids[3]
        end

        @testset "getindex(genome, Vector{String}) – no matches → empty sub-table" begin
            sub = tab[sp.genome, ["nope_1", "nope_2"]]
            @test sub isa TabularData
            @test size(sub.table) == (0, 2)
            @test isempty(sub.samples)
            @test sub.variables == tab.variables
        end

        # --- DataFrame conversion --------------------------------------------
        @testset "DataFrame(tab, genome) – columns, names and values" begin
            out = DataFrame(tab, sp.genome)
            @test names(out) == ["ID", "v1", "v2"]
            @test nrow(out) == 4
            @test out.ID[1:3] == gene_ids
            @test out.v1 == [10.0, 20.0, 30.0, 40.0]
            @test out.v2 == [11.0, 21.0, 31.0, 41.0]
        end

        @testset "DataFrame(tab, genome) – unmatched sample becomes missing" begin
            out = DataFrame(tab, sp.genome)
            @test eltype(out.ID) == Union{Missing, String}
            @test ismissing(out.ID[4])
            @test count(ismissing, out.ID) == 1
        end
    end  # TabularData indexing & DataFrame conversion

    # =========================================================================
    @testset "calculate_frequency" begin
        # Build a BedData from `scaffold => [(start, end), ...]` pairs.
        make_bed(pairs...) = begin
            scaffolds = Dict{String, IntervalMeta64}()
            for (name, ivs) in pairs
                tree = IntervalMeta64()
                for (s, e) in ivs
                    push!(tree, IntervalValue(UInt32(s), UInt32(e), UInt64(0)))
                end
                scaffolds[name] = tree
            end
            BedData(scaffolds)
        end

        @testset "per-base counts across measurements (UInt8)" begin
            m1 = make_bed("chr1" => [(1, 5), (10, 12)])
            m2 = make_bed("chr1" => [(3, 11)])
            freq = calculate_frequency([m1, m2])

            @test freq isa Dict
            @test haskey(freq, "chr1")
            v = freq["chr1"]
            @test eltype(v) == UInt8
            @test length(v) == 12
            # base:      1  2  3  4  5  6  7  8  9 10 11 12
            @test Vector(v) == UInt8[1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1]
        end

        @testset "merge default merges within-measurement overlaps" begin
            m = make_bed("chr1" => [(1, 5), (3, 8)])

            merged = calculate_frequency([m])                # merge = true (default)
            @test Vector(merged["chr1"]) == UInt8[1, 1, 1, 1, 1, 1, 1, 1]
            @test maximum(merged["chr1"]) == 1

            unmerged = calculate_frequency([m]; merge = false)
            @test Vector(unmerged["chr1"]) == UInt8[1, 1, 2, 2, 2, 1, 1, 1]
        end

        @testset "multiple scaffolds with independent coverage" begin
            m1 = make_bed("chrA" => [(1, 3)], "chrB" => [(5, 6)])
            m2 = make_bed("chrA" => [(2, 4)])
            freq = calculate_frequency([m1, m2])

            @test Set(keys(freq)) == Set(["chrA", "chrB"])
            @test Vector(freq["chrA"]) == UInt8[1, 2, 2, 1]        # length 4
            @test Vector(freq["chrB"]) == UInt8[0, 0, 0, 0, 1, 1]  # length 6
        end

        @testset "element type widens to UInt16 past 255 measurements" begin
            measurements = [make_bed("chr1" => [(1, 2)]) for _ in 1:256]
            freq = calculate_frequency(measurements)
            @test eltype(freq["chr1"]) == UInt16
            @test Vector(freq["chr1"]) == UInt16[256, 256]
        end
    end  # calculate_frequency
end
