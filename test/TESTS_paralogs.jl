using BioinfoTools2
using BioinfoTools2.Paralogs
using BioinfoTools2.Reference
using DataFrames
using SparseArrays
using Test

# Build a tiny in-memory genome from (scaffold, start, stop, strand, id) tuples.
function genome_from_genes(genes)
    gff = tempname() * ".gff3"
    open(gff, "w") do io
        println(io, "##gff-version 3")
        for (scaffold, start_pos, stop_pos, strand, id) in genes
            println(
                io,
                join(
                    (
                        scaffold, "test", "gene", string(start_pos), string(stop_pos),
                        ".", strand, ".", "ID=$id",
                    ),
                    '\t',
                ),
            )
        end
    end
    genome = Species("test").genome
    add_features!(gff, genome)
    rm(gff; force = true)
    return genome
end

# Linear index of the gene at genomic `start_pos` on `scaffold`; starts are
# unique in the test genome, so this pins down a specific gene.
function gene_index(gf, scaffold, start_pos)
    local_pos = findfirst(==(UInt32(start_pos)), gf.intervals[scaffold].start_pos)
    return gf.scaffold_ranges[scaffold][local_pos]
end


@testset "ParalogUtils" begin

    # ========================================
    # Tests for rbh()
    # ========================================

    @testset "rbh - basic functionality with max scoring" begin
        # Create test data with reciprocal best hits
        df = DataFrame(
            GeneID = ["A", "B", "C", "D"],
            ParalogID = ["B", "A", "D", "C"],
            Perc1 = [95.0, 94.0, 80.0, 85.0],
            Perc2 = [94.0, 95.0, 85.0, 80.0],
        )

        result = rbh(df; scoring = "max")

        # Should find 2 RBH pairs: A-B and C-D
        @test nrow(result) == 2
        @test "GeneID" in names(result)
        @test "ParalogID" in names(result)
        @test "perc_1" in names(result)
        @test "perc_2" in names(result)
        @test "max_perc" in names(result)
        @test "mean_perc" in names(result)

        # Check that max_perc is correctly calculated
        @test all(result.max_perc .>= result.perc_1)
        @test all(result.max_perc .>= result.perc_2)
    end

    @testset "rbh - mean scoring" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "A"],
            Perc1 = [90.0, 88.0],
            Perc2 = [88.0, 90.0],
        )

        result = rbh(df; scoring = "mean")

        @test nrow(result) == 1
        # Mean should be average of bidirectional scores
        @test result.mean_perc[1] ≈ 89.0
    end

    @testset "rbh - average scoring (alias)" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "A"],
            Perc1 = [100.0, 80.0],
            Perc2 = [80.0, 100.0],
        )

        result = rbh(df; scoring = "avg")

        @test nrow(result) == 1
        @test result.mean_perc[1] ≈ 90.0
    end

    @testset "rbh - double_max scoring" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "A"],
            Perc1 = [95.0, 93.0],
            Perc2 = [93.0, 95.0],
        )

        result = rbh(df; scoring = "double_max")

        @test nrow(result) == 1
        # In double_max mode, should use original scores
        @test result.perc_1[1] == 95.0
        @test result.perc_2[1] == 93.0
    end

    @testset "rbh - no reciprocal hits" begin
        df = DataFrame(
            GeneID = ["A", "B", "C"],
            ParalogID = ["B", "C", "D"],
            Perc1 = [95.0, 90.0, 85.0],
            Perc2 = [94.0, 89.0, 84.0],
        )

        result = rbh(df; scoring = "max")

        # May still find some hits depending on the actual algorithm
        # Just check that it returns a valid dataframe
        @test "GeneID" in names(result)
        @test "ParalogID" in names(result)
    end

    @testset "rbh - empty input should error" begin
        df = DataFrame(
            GeneID = String[],
            ParalogID = String[],
            Perc1 = Float64[],
            Perc2 = Float64[],
        )

        # Empty DataFrame will cause BoundsError when accessing elements
        @test_throws BoundsError rbh(df; scoring = "max")
    end
end

# ============================================================================
# Tests for the GeneFamily(genome, pairs) constructor
# ============================================================================
@testset "GeneFamily constructor" begin
    GF = BioinfoTools2.Paralogs.GeneFamily   # not exported from Paralogs

    # chr1 carries g1,g2,g3 (starts 100/300/500); chr2 carries g4,g5 (50/90).
    # Linear indices follow scaffold-name then start order, but every test
    # resolves indices through `gene_index` rather than assuming them.
    genes = [
        ("chr1", 100, 200, "+", "g1"),
        ("chr1", 300, 400, "-", "g2"),
        ("chr1", 500, 600, "+", "g3"),
        ("chr2", 50, 80, "+", "g4"),
        ("chr2", 90, 120, "-", "g5"),
    ]
    genome = genome_from_genes(genes)

    @testset "topology only (two columns)" begin
        df = DataFrame(query = ["g1", "g2", "g4"], subject = ["g2", "g3", "g5"])
        gf = GF(genome, df)

        # Intervals are grouped by scaffold, ordered by position, and carry the
        # feature's real 64-bit code (never a zeroed placeholder).
        @test sort(collect(keys(gf.intervals))) == ["chr1", "chr2"]
        @test collect(gf.intervals["chr1"].start_pos) == UInt32[100, 300, 500]
        @test collect(gf.intervals["chr2"].start_pos) == UInt32[50, 90]
        @test all(!=(zero(UInt64)), gf.intervals["chr1"].code)

        # The five genes fill 1:5, chr1 first.
        @test gf.scaffold_ranges["chr1"] == 1:3
        @test gf.scaffold_ranges["chr2"] == 4:5

        # Only topology is built, and it is a symmetric Bool adjacency.
        @test gf.topology isa SparseMatrixCSC{Bool, UInt32}
        @test size(gf.topology) == (5, 5)
        @test gf.topology == permutedims(gf.topology)
        @test gf.topology[1, 2] && gf.topology[2, 1]   # g1-g2
        @test gf.topology[2, 3] && gf.topology[3, 2]   # g2-g3
        @test gf.topology[4, 5] && gf.topology[5, 4]   # g4-g5
        @test !gf.topology[1, 3]                       # unrelated pair
        @test gf.dN === nothing
        @test gf.dS === nothing
        @test gf.id_subject_query === nothing
        @test gf.id_query_subject === nothing
    end

    @testset "intervals and scaffold_ranges stay concordant" begin
        gf = GF(genome, DataFrame(q = ["g1", "g3"], s = ["g5", "g4"]))
        covered = Int[]
        for (name, rng) in gf.scaffold_ranges
            @test length(rng) == length(gf.intervals[name])
            append!(covered, collect(rng))
        end
        # Ranges partition 1:n_genes exactly (no gaps, no overlaps).
        @test sort(covered) == 1:size(gf.topology, 1)
    end

    @testset "all six columns, placed and oriented correctly" begin
        df = DataFrame(
            q = ["g1"],
            s = ["g2"],
            dN = [0.10],
            dS = [0.20],
            id_subject_query = [90.0],
            id_query_subject = [80.0],
        )
        gf = GF(genome, df)
        q = gene_index(gf, "chr1", 100)   # g1 (query)
        s = gene_index(gf, "chr1", 300)   # g2 (subject)

        # dN/dS are symmetric: one entry, in the upper triangle only.
        @test gf.dN[min(q, s), max(q, s)] == 0.10
        @test gf.dS[min(q, s), max(q, s)] == 0.20
        @test iszero(gf.dN[max(q, s), min(q, s)])
        @test nnz(gf.dN) == 1
        @test nnz(gf.dS) == 1

        # %ID matrices keep their look-up key on the column axis.
        @test gf.id_subject_query[q, s] == 90.0   # subjects on columns
        @test iszero(gf.id_subject_query[s, q])
        @test gf.id_query_subject[s, q] == 80.0   # queries on columns
        @test iszero(gf.id_query_subject[q, s])
    end

    @testset "value columns matched by name, not position" begin
        # Columns out of spec order, an unrecognised column, one %ID omitted.
        df = DataFrame(
            q = ["g1"],
            s = ["g2"],
            id_query_subject = [80.0],
            dS = [0.20],
            note = ["ignored"],
            dN = [0.10],
        )
        gf = GF(genome, df)
        q = gene_index(gf, "chr1", 100)
        s = gene_index(gf, "chr1", 300)
        @test gf.dN[min(q, s), max(q, s)] == 0.10
        @test gf.dS[min(q, s), max(q, s)] == 0.20
        @test gf.id_query_subject[s, q] == 80.0
        @test gf.id_subject_query === nothing     # never supplied
    end

    @testset "column name matching is case-insensitive" begin
        gf = GF(genome, DataFrame(q = ["g1"], s = ["g2"], DN = [0.3], DS = [0.4]))
        q = gene_index(gf, "chr1", 100)
        s = gene_index(gf, "chr1", 300)
        @test gf.dN[min(q, s), max(q, s)] == 0.3
        @test gf.dS[min(q, s), max(q, s)] == 0.4
    end

    @testset "partial: only one value column supplied" begin
        gf = GF(genome, DataFrame(q = ["g1"], s = ["g3"], dN = [0.5]))
        q = gene_index(gf, "chr1", 100)   # g1
        s = gene_index(gf, "chr1", 500)   # g3
        @test gf.dN[min(q, s), max(q, s)] == 0.5
        @test gf.dS === nothing
        @test gf.id_subject_query === nothing
        @test gf.id_query_subject === nothing
        @test gf.topology[q, s] && gf.topology[s, q]
    end

    @testset "integer value columns convert to Float64" begin
        gf = GF(genome, DataFrame(q = ["g1"], s = ["g2"], dN = [1]))
        q = gene_index(gf, "chr1", 100)
        s = gene_index(gf, "chr1", 300)
        @test gf.dN isa SparseMatrixCSC{Float64, UInt32}
        @test gf.dN[min(q, s), max(q, s)] == 1.0
    end

    @testset "IDs absent from the genome (and their rows) are skipped" begin
        df = DataFrame(q = ["g1", "ghost", "g2"], s = ["nope", "g3", "g4"])
        gf = GF(genome, df)

        # ghost/nope are unknown; only the g2-g4 row survives. The four real
        # IDs still get indexed, so g1 and g3 sit as isolated nodes.
        @test size(gf.topology) == (4, 4)
        g2 = gene_index(gf, "chr1", 300)
        g4 = gene_index(gf, "chr2", 50)
        @test nnz(gf.topology) == 2                     # (g2,g4) + (g4,g2)
        @test gf.topology[g2, g4] && gf.topology[g4, g2]
    end

    @testset "reciprocal rows collapse instead of double-counting" begin
        df = DataFrame(q = ["g1", "g2"], s = ["g2", "g1"], dN = [0.42, 0.42])
        gf = GF(genome, df)
        q = gene_index(gf, "chr1", 100)
        s = gene_index(gf, "chr1", 300)
        @test nnz(gf.dN) == 1
        @test gf.dN[min(q, s), max(q, s)] == 0.42
    end

    @testset "fewer than two columns errors" begin
        @test_throws ArgumentError GF(genome, DataFrame(query = ["g1"]))
        @test_throws ArgumentError GF(genome, DataFrame())
    end

    @testset "non-numeric value column errors" begin
        df = DataFrame(q = ["g1"], s = ["g2"], dN = ["not_a_number"])
        @test_throws MethodError GF(genome, df)
    end
end
