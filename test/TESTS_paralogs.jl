using BioinfoTools2
using BioinfoTools2.Homologs.Paralogs
using BioinfoTools2.Reference
using DataFrames
using Graphs
using SimpleWeightedGraphs
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
                        scaffold,
                        "test",
                        "gene",
                        string(start_pos),
                        string(stop_pos),
                        ".",
                        strand,
                        ".",
                        "ID=$id",
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
function gene_index(pg, scaffold, start_pos)
    local_pos = findfirst(==(UInt32(start_pos)), pg.intervals[scaffold].start_pos)
    return pg.scaffold_ranges[scaffold][local_pos]
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

    @testset "rbh - empty input yields an empty (but shaped) result" begin
        df = DataFrame(
            GeneID = String[],
            ParalogID = String[],
            Perc1 = Float64[],
            Perc2 = Float64[],
        )

        result = rbh(df; scoring = "max")
        @test nrow(result) == 0
        @test names(result) ==
              ["GeneID", "ParalogID", "perc_1", "perc_2", "max_perc", "mean_perc"]
    end

    @testset "rbh - column type and arity validation" begin
        @test_throws ArgumentError rbh(
            DataFrame(GeneID = ["A"], ParalogID = ["B"], Perc1 = [95.0]),
        )
        @test_throws ArgumentError rbh(
            DataFrame(GeneID = [1], ParalogID = ["B"], Perc1 = [95.0], Perc2 = [94.0]),
        )
        @test_throws ArgumentError rbh(
            DataFrame(GeneID = ["A"], ParalogID = ["B"], Perc1 = ["x"], Perc2 = [94.0]),
        )
        @test_throws ArgumentError rbh(
            DataFrame(GeneID = ["A"], ParalogID = ["B"], Perc1 = [95.0], Perc2 = [94.0]);
            scoring = "not_a_scoring",
        )
    end

    @testset "rbh_ds - lowest dS wins, unscored pairs ignored" begin
        # a-b are each other's closest pair (dS 0.05); c-d likewise (0.30).
        # a-d is a distant pair that must not beat either.
        df = DataFrame(
            GeneID = ["a", "c", "a"],
            ParalogID = ["b", "d", "d"],
            dS = [0.05, 0.30, 0.90],
        )

        result = rbh_ds(df)
        @test names(result) == ["GeneID", "ParalogID", "ds", "min_ds"]
        @test nrow(result) == 2

        pairs = Set(Set([row.GeneID, row.ParalogID]) for row in eachrow(result))
        @test pairs == Set([Set(["a", "b"]), Set(["c", "d"])])
        @test sort(result.ds) == [0.05, 0.30]
        @test result.ds == result.min_ds

        # `scoring = "ds"` routes to the same result.
        @test rbh(df; scoring = "ds") == result
    end
end

# ============================================================================
# Tests for the ParalogGroup(genome, pairs) constructor
# ============================================================================
@testset "ParalogGroup constructor" begin
    PG = BioinfoTools2.Homologs.Paralogs.ParalogGroup   # not exported from Paralogs

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
        pg = PG(genome, df)

        # Intervals are grouped by scaffold, ordered by position, and carry the
        # feature's real 64-bit code (never a zeroed placeholder).
        @test sort(collect(keys(pg.intervals))) == ["chr1", "chr2"]
        @test collect(pg.intervals["chr1"].start_pos) == UInt32[100, 300, 500]
        @test collect(pg.intervals["chr2"].start_pos) == UInt32[50, 90]
        @test all(!=(zero(UInt64)), pg.intervals["chr1"].code)

        # The five genes fill 1:5, chr1 first.
        @test pg.scaffold_ranges["chr1"] == 1:3
        @test pg.scaffold_ranges["chr2"] == 4:5

        # Only topology is built, and it is a symmetric Bool adjacency.
        @test pg.topology isa SparseMatrixCSC{Bool,UInt32}
        @test size(pg.topology) == (5, 5)
        @test pg.topology == permutedims(pg.topology)
        @test pg.topology[1, 2] && pg.topology[2, 1]   # g1-g2
        @test pg.topology[2, 3] && pg.topology[3, 2]   # g2-g3
        @test pg.topology[4, 5] && pg.topology[5, 4]   # g4-g5
        @test !pg.topology[1, 3]                       # unrelated pair
        @test pg.dN === nothing
        @test pg.dS === nothing
        @test pg.id_subject_query === nothing
        @test pg.id_query_subject === nothing
    end

    @testset "intervals and scaffold_ranges stay concordant" begin
        pg = PG(genome, DataFrame(q = ["g1", "g3"], s = ["g5", "g4"]))
        covered = Int[]
        lengths_match = Bool[]
        for (name, rng) in pg.scaffold_ranges
            push!(lengths_match, length(rng) == length(pg.intervals[name]))
            append!(covered, collect(rng))
        end
        @test all(lengths_match)
        # Ranges partition 1:n_genes exactly (no gaps, no overlaps).
        @test sort(covered) == 1:size(pg.topology, 1)
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
        pg = PG(genome, df)
        q = gene_index(pg, "chr1", 100)   # g1 (query)
        s = gene_index(pg, "chr1", 300)   # g2 (subject)

        # dN/dS are symmetric: one entry, in the upper triangle only.
        @test pg.dN[min(q, s), max(q, s)] == 0.10
        @test pg.dS[min(q, s), max(q, s)] == 0.20
        @test iszero(pg.dN[max(q, s), min(q, s)])
        @test nnz(pg.dN) == 1
        @test nnz(pg.dS) == 1

        # %ID matrices keep their look-up key on the column axis.
        @test pg.id_subject_query[q, s] == 90.0   # subjects on columns
        @test iszero(pg.id_subject_query[s, q])
        @test pg.id_query_subject[s, q] == 80.0   # queries on columns
        @test iszero(pg.id_query_subject[q, s])
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
        pg = PG(genome, df)
        q = gene_index(pg, "chr1", 100)
        s = gene_index(pg, "chr1", 300)
        @test pg.dN[min(q, s), max(q, s)] == 0.10
        @test pg.dS[min(q, s), max(q, s)] == 0.20
        @test pg.id_query_subject[s, q] == 80.0
        @test pg.id_subject_query === nothing     # never supplied
    end

    @testset "column name matching is case-insensitive" begin
        pg = PG(genome, DataFrame(q = ["g1"], s = ["g2"], DN = [0.3], DS = [0.4]))
        q = gene_index(pg, "chr1", 100)
        s = gene_index(pg, "chr1", 300)
        @test pg.dN[min(q, s), max(q, s)] == 0.3
        @test pg.dS[min(q, s), max(q, s)] == 0.4
    end

    @testset "partial: only one value column supplied" begin
        pg = PG(genome, DataFrame(q = ["g1"], s = ["g3"], dN = [0.5]))
        q = gene_index(pg, "chr1", 100)   # g1
        s = gene_index(pg, "chr1", 500)   # g3
        @test pg.dN[min(q, s), max(q, s)] == 0.5
        @test pg.dS === nothing
        @test pg.id_subject_query === nothing
        @test pg.id_query_subject === nothing
        @test pg.topology[q, s] && pg.topology[s, q]
    end

    @testset "integer value columns convert to Float64" begin
        pg = PG(genome, DataFrame(q = ["g1"], s = ["g2"], dN = [1]))
        q = gene_index(pg, "chr1", 100)
        s = gene_index(pg, "chr1", 300)
        @test pg.dN isa SparseMatrixCSC{Float64,UInt32}
        @test pg.dN[min(q, s), max(q, s)] == 1.0
    end

    @testset "IDs absent from the genome (and their rows) are skipped" begin
        df = DataFrame(q = ["g1", "ghost", "g2"], s = ["nope", "g3", "g4"])
        pg = PG(genome, df)

        # ghost/nope are unknown, so only the g2-g4 row survives. g1 and g3 lose
        # their only partner with it, and never enter the group.
        @test sort(collect(keys(pg.id_to_index))) == ["g2", "g4"]
        @test size(pg.topology) == (2, 2)
        g2 = gene_index(pg, "chr1", 300)
        g4 = gene_index(pg, "chr2", 50)
        @test nnz(pg.topology) == 2                     # (g2,g4) + (g4,g2)
        @test pg.topology[g2, g4] && pg.topology[g4, g2]
    end

    @testset "rows with 0% identity in every direction are dropped" begin
        df = DataFrame(
            q = ["g1", "g2"],
            s = ["g2", "g3"],
            id_subject_query = [90.0, 0.0],
            id_query_subject = [80.0, 0.0],
        )
        local pg
        pg = @test_logs (:warn, r"dropped 1 row with 0% identity") PG(genome, df)
        # g3 was only ever paired through the dropped row.
        @test sort(collect(keys(pg.id_to_index))) == ["g1", "g2"]
        @test size(pg.topology) == (2, 2)

        # One non-zero direction is enough to keep the pair.
        one_way = PG(
            genome,
            DataFrame(
                q = ["g1"],
                s = ["g2"],
                id_subject_query = [0.0],
                id_query_subject = [5.0],
            ),
        )
        @test size(one_way.topology) == (2, 2)

        # With a single %ID column supplied, that column alone decides.
        local empty_group
        empty_group = @test_logs (:warn, r"0% identity") PG(
            genome,
            DataFrame(q = ["g1"], s = ["g2"], id_query_subject = [0.0]),
        )
        @test size(empty_group.topology) == (0, 0)
        @test isempty(empty_group.id_to_index)

        # No %ID column at all: the pair is taken on trust, dS 0.0 included.
        trusted = @test_logs PG(genome, DataFrame(q = ["g1"], s = ["g2"], dS = [0.0]))
        @test size(trusted.topology) == (2, 2)
    end

    @testset "reciprocal rows collapse instead of double-counting" begin
        df = DataFrame(q = ["g1", "g2"], s = ["g2", "g1"], dN = [0.42, 0.42])
        pg = PG(genome, df)
        q = gene_index(pg, "chr1", 100)
        s = gene_index(pg, "chr1", 300)
        @test nnz(pg.dN) == 1
        @test pg.dN[min(q, s), max(q, s)] == 0.42
    end

    @testset "a zero relation value is kept, and stays a real zero" begin
        # `sparse` stores an explicit 0, so a dN/dS of 0 survives construction,
        # slicing and `rbh` intact — only the weighted-graph views drop it.
        pg = PG(genome, DataFrame(q = ["g1"], s = ["g2"], dS = [0.0], dN = [0.0]))
        q = gene_index(pg, "chr1", 100)
        s = gene_index(pg, "chr1", 300)
        @test pg.topology[q, s] && pg.topology[s, q]
        @test nnz(pg.dS) == 1 && nnz(pg.dN) == 1
        @test pg.dS[min(q, s), max(q, s)] == 0.0
        @test last(findnz(pg.dS)) == [0.0]
        @test nnz(pg[["g1", "g2"]].dS) == 1     # sub-group slicing keeps it

        # A dS of 0 is the best possible rank, and round-trips as a real 0.
        result = rbh(pg)
        @test nrow(result) == 1
        @test result.dS[1] == 0.0
        @test nnz(PG(genome, result).topology) == 2

        # `m + permutedims(m)` prunes it, and a weighted graph could not hold a
        # zero-weight edge anyway.
        @test ne(dS_graph(pg)) == 0
    end

    @testset "NaN rows are dropped up front, with a warning" begin
        # A NaN in *either* matched column drops the whole row, not just that
        # field's contribution.
        df = DataFrame(
            q = ["g1", "g2", "g1"],
            s = ["g2", "g3", "g3"],
            dN = [0.1, NaN, 0.3],
            dS = [0.4, 0.5, NaN],
        )
        local pg
        pg = @test_logs (:warn, r"dropped 2 rows") PG(genome, df)

        # Only the g1-g2 row (no NaNs) survives.
        @test size(pg.topology) == (2, 2)
        @test nnz(pg.dN) == 1
        @test nnz(pg.dS) == 1
        @test !any(isnan, nonzeros(pg.dN))
        @test !any(isnan, nonzeros(pg.dS))

        # No NaNs -> no warning, nothing dropped.
        @test_logs PG(genome, DataFrame(q = ["g1"], s = ["g2"], dN = [0.1]))

        # NaN in an unrecognised (non-value) column never triggers dropping.
        pg2 = @test_logs PG(genome, DataFrame(q = ["g1"], s = ["g2"], note = [NaN]))
        @test size(pg2.topology) == (2, 2)
        @test pg2.dN === nothing
    end

    @testset "fewer than two columns errors" begin
        @test_throws ArgumentError PG(genome, DataFrame(query = ["g1"]))
        @test_throws ArgumentError PG(genome, DataFrame())
    end

    @testset "column type validation" begin
        # Query ID column holding non-strings.
        @test_throws ArgumentError PG(genome, DataFrame(q = [1, 2], s = ["g2", "g3"]))

        # Subject ID column holding non-strings.
        @test_throws ArgumentError PG(genome, DataFrame(q = ["g1", "g2"], s = [1, 2]))

        # A recognised value column holding non-numeric data.
        df_bad_dn = DataFrame(q = ["g1"], s = ["g2"], dN = ["not_a_number"])
        @test_throws ArgumentError PG(genome, df_bad_dn)

        df_bad_idsq = DataFrame(q = ["g1"], s = ["g2"], id_subject_query = ["high"])
        @test_throws ArgumentError PG(genome, df_bad_idsq)

        # The `lca` column is categorical, so it must hold strings, not numbers.
        @test_throws ArgumentError PG(
            genome,
            DataFrame(q = ["g1"], s = ["g2"], lca = [3.5]),
        )

        # Error messages name the offending column and its element type.
        try
            PG(genome, DataFrame(q = [1], s = ["g2"]))
            @test false   # unreachable
        catch err
            @test err isa ArgumentError
            @test occursin("Query ID column", err.msg)
        end
        try
            PG(genome, DataFrame(q = ["g1"], s = ["g2"], dS = ["x"]))
            @test false   # unreachable
        catch err
            @test err isa ArgumentError
            @test occursin("dS", err.msg)
        end

        # An unrecognised, non-value column is left alone regardless of type.
        pg = PG(genome, DataFrame(q = ["g1"], s = ["g2"], note = [1, 2, 3][1:1]))
        @test pg.dN === nothing
    end
end

# ============================================================================
# Tests for the `lca` relation
# ============================================================================
@testset "ParalogGroup lca" begin
    PG = BioinfoTools2.Homologs.Paralogs.ParalogGroup
    genes = [
        ("chr1", 100, 200, "+", "g1"),
        ("chr1", 300, 400, "-", "g2"),
        ("chr1", 500, 600, "+", "g3"),
    ]
    genome = genome_from_genes(genes)
    df = DataFrame(
        q = ["g1", "g2"],
        s = ["g2", "g3"],
        dS = [0.1, 0.2],
        lca = ["Hominidae", "Primates"],
    )
    pg = PG(genome, df)
    g1 = gene_index(pg, "chr1", 100)
    g2 = gene_index(pg, "chr1", 300)
    g3 = gene_index(pg, "chr1", 500)

    @testset "interned codes, stored in the upper triangle" begin
        @test pg.lca isa SparseMatrixCSC{UInt32,UInt32}
        @test sort(pg.lca_labels) == ["Hominidae", "Primates"]
        @test nnz(pg.lca) == 2
        @test pg.lca_labels[pg.lca[min(g1, g2), max(g1, g2)]] == "Hominidae"
        @test iszero(pg.lca[max(g1, g2), min(g1, g2)])   # mirror stays empty
    end

    @testset "lca_label reads either orientation, by index or by ID" begin
        @test lca_label(pg, g1, g2) == "Hominidae"
        @test lca_label(pg, g2, g1) == "Hominidae"
        @test lca_label(pg, "g2", "g3") == "Primates"
        @test lca_label(pg, g1, g3) === nothing   # never paired in the table
        @test_throws ArgumentError lca_label(pg, "g1", "ghost")
        @test_throws ArgumentError lca_label(pg, 0, g1)

        without_lca = PG(genome, DataFrame(q = ["g1"], s = ["g2"]))
        @test without_lca.lca === nothing
        @test isempty(without_lca.lca_labels)
        @test lca_label(without_lca, "g1", "g2") === nothing
    end

    @testset "repeated labels intern once; an empty label means unknown" begin
        shared = PG(
            genome,
            DataFrame(q = ["g1", "g2"], s = ["g2", "g3"], lca = ["Primates", "Primates"]),
        )
        @test shared.lca_labels == ["Primates"]
        @test nnz(shared.lca) == 2

        blank = PG(
            genome,
            DataFrame(q = ["g1", "g2"], s = ["g2", "g3"], lca = ["", "Primates"]),
        )
        @test blank.lca_labels == ["Primates"]
        @test lca_label(blank, "g1", "g2") === nothing
        # An unknown ancestor is not a missing pair: topology keeps the edge.
        @test blank.topology[blank.id_to_index["g1"], blank.id_to_index["g2"]]
    end

    @testset "surfaced by getindex, show, sub-groups and rbh" begin
        col = pg[g2]
        @test col[:lca] isa Vector{String}
        @test length(col[:lca]) == size(pg.topology, 1)
        @test col[:lca][g1] == "Hominidae"
        @test col[:lca][g3] == ""     # g2-g3 lives in column g3, not here

        @test occursin("lca", sprint(show, pg))

        sub = pg[["g1", "g2"]]
        @test lca_label(sub, "g1", "g2") == "Hominidae"

        # rbh keeps the winning pair's ancestor, so its output round-trips.
        result = rbh(pg)
        @test "lca" in names(result)
        @test result.lca[1] == "Hominidae"   # g1-g2 wins on the lower dS
        round_trip = PG(genome, result)
        @test lca_label(round_trip, result.query[1], result.subject[1]) == result.lca[1]
    end
end

# ============================================================================
# Tests for the `lca_depth` relation
# ============================================================================
@testset "ParalogGroup lca_depth" begin
    PG = BioinfoTools2.Homologs.Paralogs.ParalogGroup
    genes = [
        ("chr1", 100, 200, "+", "g1"),
        ("chr1", 300, 400, "-", "g2"),
        ("chr1", 500, 600, "+", "g3"),
    ]
    genome = genome_from_genes(genes)
    df = DataFrame(
        q = ["g1", "g1"],
        s = ["g2", "g3"],
        dS = [0.5, 0.1],
        lca = ["Hominidae", "Primates"],
        lca_depth = [2.0, 1.0],
    )
    pg = PG(genome, df)
    g1 = gene_index(pg, "chr1", 100)
    g2 = gene_index(pg, "chr1", 300)
    g3 = gene_index(pg, "chr1", 500)

    @testset "numeric and symmetric, stored once in the upper triangle" begin
        @test pg.lca_depth isa SparseMatrixCSC{Float64,UInt32}
        @test nnz(pg.lca_depth) == 2
        @test pg.lca_depth[min(g1, g2), max(g1, g2)] == 2.0
        @test iszero(pg.lca_depth[max(g1, g2), min(g1, g2)])   # mirror stays empty

        # Integer depths convert, like any other numeric relation column.
        ints = PG(genome, DataFrame(q = ["g1"], s = ["g2"], lca_depth = [3]))
        iq, is = gene_index(ints, "chr1", 100), gene_index(ints, "chr1", 300)
        @test ints.lca_depth[min(iq, is), max(iq, is)] == 3.0

        # A NaN depth drops its row, like any other numeric relation column.
        local dropped
        dropped = @test_logs (:warn, r"dropped 1 row") PG(
            genome,
            DataFrame(q = ["g1", "g1"], s = ["g2", "g3"], lca_depth = [NaN, 1.0]),
        )
        @test nnz(dropped.lca_depth) == 1
    end

    @testset "orders a gene's partners, and pairs against each other" begin
        # Row plus column is a gene's complete symmetric view; the diagonal is
        # never set, so nothing is double-counted.
        depths = Vector(pg.lca_depth[:, g1]) .+ Vector(pg.lca_depth[g1, :])
        @test depths[g2] == 2.0
        @test depths[g3] == 1.0

        # Deeper is more recent, so ascending depth is oldest duplication first.
        partners = [g2, g3]
        @test partners[sortperm(depths[partners])] == [g3, g2]
    end

    @testset "surfaced by getindex, show, sub-groups and rbh" begin
        col = pg[g2]
        @test haskey(col, :lca_depth)
        @test col[:lca_depth][g1] == 2.0
        @test occursin("lca_depth", sprint(show, pg))

        sub = pg[["g1", "g2"]]
        sg1, sg2 = gene_index(sub, "chr1", 100), gene_index(sub, "chr1", 300)
        @test sub.lca_depth[min(sg1, sg2), max(sg1, sg2)] == 2.0

        # g1-g3 wins the component on the lower dS, and carries its depth out.
        result = rbh(pg)
        @test "lca_depth" in names(result)
        @test result.lca_depth[1] == 1.0
        round_trip = PG(genome, result)
        @test nnz(round_trip.lca_depth) == 1
    end

    @testset "never ranks rbh unless asked for" begin
        depth_only = PG(genome, DataFrame(q = ["g1"], s = ["g2"], lca_depth = [2.0]))
        @test_throws ArgumentError rbh(depth_only)
        @test nrow(rbh(depth_only; levels = [:lca_depth])) == 1
    end
end

# ============================================================================
# Tests for add_duplications_of (OrthoFinder Duplications.tsv ingestion)
# ============================================================================
@testset "add_duplications_of" begin
    PG = BioinfoTools2.Homologs.Paralogs.ParalogGroup
    depths_of = BioinfoTools2.Homologs.Paralogs._newick_depths

    # Shaped like OrthoFinder's own tree. For `human`'s ancestors the depths are
    # N0 = 1, N1 = 2, N2 = 3, and the terminal `human` = 4.
    tree = "(mouse:0.08,(macaque:0.02,(human:0.01,(bonobo:0.004,chimp:0.004)N3:0.003)N2:0.005)N1:0.08)N0;"

    # An OrthoFinder-shaped results directory, returning its Duplications.tsv.
    function duplications_file(rows; newick = tree)
        root = mktempdir()
        mkpath(joinpath(root, "Species_Tree"))
        mkpath(joinpath(root, "Gene_Duplication_Events"))
        write(joinpath(root, "Species_Tree", "SpeciesTree_rooted_node_labels.txt"), newick)
        path = joinpath(root, "Gene_Duplication_Events", "Duplications.tsv")
        open(path, "w") do io
            println(
                io,
                join(
                    (
                        "Orthogroup",
                        "Species Tree Node",
                        "Gene Tree Node",
                        "Support",
                        "Type",
                        "Genes 1",
                        "Genes 2",
                    ),
                    '\t',
                ),
            )
            for row in rows
                println(io, join(row, '\t'))
            end
        end
        return path
    end

    genes = [
        ("chr1", 100, 200, "+", "g1"),
        ("chr1", 300, 400, "-", "g2"),
        ("chr1", 500, 600, "+", "g3"),
        ("chr1", 700, 800, "+", "g4"),
    ]
    genome = genome_from_genes(genes)
    # g1-g2, g1-g3 and g3-g4 are pairs; only two of them are duplications the
    # file records.
    pairs = DataFrame(
        q = ["g1", "g1", "g3"],
        s = ["g2", "g3", "g4"],
        dS = [0.5, 0.7, 0.9],
        id_subject_query = [90.0, 80.0, 70.0],
        id_query_subject = [91.0, 81.0, 71.0],
    )
    pg = PG(genome, pairs)

    rows = [
        # A non-terminal node, with a non-target species mixed into both sides.
        ("OG0", "N2", "n1", "1.0", "Non-Terminal", "human_g1, mouse_g9", "human_g2"),
        # A terminal node is labelled with the species name itself.
        ("OG0", "human", "n2", "1.0", "Terminal", "human_g3", "human_g4"),
        # Another species' duplication, which must not contribute anything.
        ("OG1", "N1", "n3", "1.0", "Non-Terminal", "mouse_g1", "mouse_g2"),
    ]

    @testset "depths come from the tree, never the node number" begin
        depths = depths_of(tree)
        @test depths["N0"] == 1
        @test depths["N1"] == 2
        @test depths["N2"] == 3
        @test depths["human"] == 4
        @test depths["bonobo"] == depths["chimp"] == 5
        @test depths["mouse"] == 2

        # OrthoFinder numbers nodes by traversal order, so a number can run
        # against depth; the nesting is what counts.
        reversed = depths_of("(a,(b,(c,d)N1)N2)N3;")
        @test reversed["N3"] == 1 && reversed["N2"] == 2 && reversed["N1"] == 3
    end

    @testset "annotates matched pairs and drops the rest" begin
        local result
        result = @test_logs (:warn, r"dropped 1 pair") add_duplications_of(
            pg,
            duplications_file(rows),
            "human",
        )

        # Only the two recorded duplications remain, all four genes with them.
        @test nnz(result.topology) == 4
        @test sort(collect(keys(result.id_to_index))) == ["g1", "g2", "g3", "g4"]
        @test lca_label(result, "g1", "g2") == "N2"
        @test lca_label(result, "g3", "g4") == "human"
        @test sort(result.lca_labels) == ["N2", "human"]

        depth(a, b) =
            result.lca_depth[minmax(result.id_to_index[a], result.id_to_index[b])...]
        @test depth("g1", "g2") == 3.0
        @test depth("g3", "g4") == 4.0

        # The unrecorded pair loses its edge *and* the values it carried.
        @test !result.topology[result.id_to_index["g1"], result.id_to_index["g3"]]
        @test nnz(result.dS) == 2
        @test result.dS[minmax(result.id_to_index["g1"], result.id_to_index["g2"])...] ==
              0.5
        @test result.id_subject_query[result.id_to_index["g1"], result.id_to_index["g2"]] ==
              90.0

        # Every surviving pair has a depth, so ranking by it now works.
        @test nrow(rbh(result; levels = [:lca_depth, :dS])) == 2

        # `pg` itself is untouched.
        @test nnz(pg.topology) == 6
        @test pg.lca === nothing
    end

    @testset "a gene left without a partner goes too" begin
        # Only g1-g2 is recorded, so g3 and g4 lose their only pairs.
        local result
        result = @test_logs (:warn, r"dropped 2 pair\(s\).*2 gene") add_duplications_of(
            pg,
            duplications_file(rows[1:1]),
            "human",
        )
        @test sort(collect(keys(result.id_to_index))) == ["g1", "g2"]
        @test sort(collect(keys(result.intervals))) == ["chr1"]
        @test length(result.intervals["chr1"]) == 2
    end

    @testset "overwrites any LCA the group already carried" begin
        stale = PG(
            genome,
            DataFrame(
                q = ["g1"],
                s = ["g2"],
                dS = [0.5],
                lca = ["Stale"],
                lca_depth = [99.0],
            ),
        )
        result = add_duplications_of(stale, duplications_file(rows), "human")
        @test lca_label(result, "g1", "g2") == "N2"
        @test result.lca_labels == ["N2"]     # the stale label is gone entirely
        @test result.lca_depth[minmax(
            result.id_to_index["g1"],
            result.id_to_index["g2"],
        )...] == 3.0
    end

    @testset "species prefixes are matched whole, underscores and all" begin
        # `pongo_abelii_g1` must not be read as species `pongo`.
        pongo_tree = "(mouse:0.1,(pongo:0.1,pongo_abelii:0.1)N1:0.1)N0;"
        pongo_rows = [(
            "OG0",
            "pongo_abelii",
            "n1",
            "1.0",
            "Terminal",
            "pongo_abelii_g1, pongo_g3",
            "pongo_abelii_g2",
        )]
        path = duplications_file(pongo_rows; newick = pongo_tree)

        result =
            @test_logs (:warn, r"dropped") add_duplications_of(pg, path, "pongo_abelii")
        @test sort(collect(keys(result.id_to_index))) == ["g1", "g2"]
        @test lca_label(result, "g1", "g2") == "pongo_abelii"

        # `pongo`'s own row names only g3, which pairs with nothing here.
        @test_throws ArgumentError add_duplications_of(pg, path, "pongo")
    end

    @testset "min_support and transform" begin
        weak = [("OG0", "N2", "n1", "0.3", "Non-Terminal", "human_g1", "human_g2")]
        path = duplications_file(weak)
        @test lca_label(add_duplications_of(pg, path, "human"), "g1", "g2") == "N2"
        @test_throws ArgumentError add_duplications_of(pg, path, "human"; min_support = 0.5)

        # Names needing a fix-up before they match the genome's IDs.
        renamed = [("OG0", "N2", "n1", "1.0", "Non-Terminal", "human_G1", "human_G2")]
        renamed_path = duplications_file(renamed)
        @test_throws ArgumentError add_duplications_of(pg, renamed_path, "human")
        result = add_duplications_of(pg, renamed_path, "human"; transform = lowercase)
        @test lca_label(result, "g1", "g2") == "N2"
    end

    @testset "input problems are errors, not silence" begin
        path = duplications_file(rows)

        # A species the tree doesn't hold.
        @test_throws ArgumentError add_duplications_of(pg, path, "zebrafish")

        # No species tree to take depths from.
        orphan = joinpath(mktempdir(), "Duplications.tsv")
        cp(path, orphan)
        @test_throws ArgumentError add_duplications_of(pg, orphan, "human")
        @test lca_label(
            add_duplications_of(
                pg,
                orphan,
                "human";
                species_tree = joinpath(
                    dirname(dirname(path)),
                    "Species_Tree",
                    "SpeciesTree_rooted_node_labels.txt",
                ),
            ),
            "g1",
            "g2",
        ) == "N2"

        # A file that isn't a Duplications.tsv at all.
        wrong = duplications_file(rows)
        write(wrong, "not\ta\tduplications\tfile\n")
        @test_throws ArgumentError add_duplications_of(pg, wrong, "human")
    end
end

# ============================================================================
# Tests for rbh's `levels` ranking precedence
# ============================================================================
@testset "rbh ranking levels" begin
    PG = BioinfoTools2.Homologs.Paralogs.ParalogGroup
    genes = [
        ("chr1", 100, 200, "+", "g1"),
        ("chr1", 300, 400, "-", "g2"),
        ("chr1", 500, 600, "+", "g3"),
    ]
    genome = genome_from_genes(genes)

    # One component. g1-g2 is the more recent duplication (deeper LCA) but the
    # higher dS; g1-g3 the reverse, so the two levels disagree outright.
    df = DataFrame(
        q = ["g1", "g1"],
        s = ["g2", "g3"],
        dS = [0.9, 0.1],
        lca_depth = [3.0, 2.0],
    )
    pg = PG(genome, df)
    winner(result) = Set([result.query[1], result.subject[1]])

    @testset "lca_depth ranks descending, and only when asked for" begin
        @test winner(rbh(pg)) == Set(["g1", "g3"])                      # dS decides
        @test winner(rbh(pg; levels = [:lca_depth])) == Set(["g1", "g2"])
        @test winner(rbh(pg; levels = [:lca_depth, :dS])) == Set(["g1", "g2"])
        @test winner(rbh(pg; levels = [:dS, :lca_depth])) == Set(["g1", "g3"])
    end

    @testset "a later level breaks an earlier tie, without warning" begin
        tied_depths = PG(
            genome,
            DataFrame(
                q = ["g1", "g1"],
                s = ["g2", "g3"],
                dS = [0.9, 0.1],
                lca_depth = [2.0, 2.0],
            ),
        )
        local result
        result = @test_logs rbh(tied_depths; levels = [:lca_depth, :dS])
        @test winner(result) == Set(["g1", "g3"])

        # Tying on *every* level is what warrants the warning.
        @test_logs (:warn, r"tied best-hit") rbh(tied_depths; levels = [:lca_depth])
    end

    @testset "a missing rank is an error, never a default" begin
        # A level `pg` doesn't carry at all.
        @test_throws ArgumentError rbh(pg; levels = [:dN])
        @test_throws ArgumentError rbh(pg; levels = [:identity])

        # An individual pair with no depth recorded (a structural zero).
        partial = PG(
            genome,
            DataFrame(
                q = ["g1", "g1"],
                s = ["g2", "g3"],
                dS = [0.9, 0.1],
                lca_depth = [2.0, 0.0],
            ),
        )
        @test_throws ArgumentError rbh(partial; levels = [:lca_depth])
        @test nrow(rbh(partial)) == 1   # dS alone still ranks it
    end

    @testset "levels are validated" begin
        @test_throws ArgumentError rbh(pg; levels = [:dS, :dS])
        @test_throws ArgumentError rbh(pg; levels = Symbol[])
        @test_throws ArgumentError rbh(pg; levels = [:not_a_level])

        # A group with nothing rankable by default still says so.
        bare = PG(genome, DataFrame(q = ["g1"], s = ["g2"]))
        @test_throws ArgumentError rbh(bare)
    end
end

# ============================================================================
# Tests for Base.show(::ParalogGroup)
# ============================================================================
@testset "ParalogGroup show" begin
    PG = BioinfoTools2.Homologs.Paralogs.ParalogGroup
    genes = [
        ("chr1", 100, 200, "+", "g1"),
        ("chr1", 300, 400, "-", "g2"),
        ("chr2", 50, 80, "+", "g3"),
    ]
    genome = genome_from_genes(genes)

    @testset "topology only" begin
        pg = PG(genome, DataFrame(q = ["g1"], s = ["g2"]))
        str = sprint(show, pg)
        @test occursin("ParalogGroup(", str)
        # Only g1/g2 (both on chr1) are referenced, so g3/chr2 never get indexed.
        @test occursin("2 genes", str)
        @test occursin("1 scaffold", str)
        @test occursin("1 pair", str)
        @test occursin("none", str)   # no relation matrices present
    end

    @testset "with relations" begin
        df = DataFrame(q = ["g1", "g2"], s = ["g2", "g3"], dN = [0.1, 0.2])
        pg = PG(genome, df)
        str = sprint(show, pg)
        @test occursin("2 pairs", str)
        @test occursin("dN", str)
        @test !occursin("dS", str)   # dS was never supplied
    end
end

# ============================================================================
# Tests for ParalogGroup getindex overloads
# ============================================================================
@testset "ParalogGroup getindex" begin
    PG = BioinfoTools2.Homologs.Paralogs.ParalogGroup
    genes = [
        ("chr1", 100, 200, "+", "g1"),
        ("chr1", 300, 400, "-", "g2"),
        ("chr1", 500, 600, "+", "g3"),
        ("chr2", 50, 80, "+", "g4"),
        ("chr2", 90, 120, "-", "g5"),
    ]
    genome = genome_from_genes(genes)
    df = DataFrame(
        q = ["g1", "g2", "g4"],
        s = ["g2", "g3", "g5"],
        dN = [0.1, 0.2, 0.3],
        dS = [0.4, 0.5, 0.6],
        id_subject_query = [90.0, 91.0, 92.0],
        id_query_subject = [80.0, 81.0, 82.0],
    )
    pg = PG(genome, df)
    g1 = gene_index(pg, "chr1", 100)
    g2 = gene_index(pg, "chr1", 300)

    @testset "pg[i::Integer] - column dict" begin
        col = pg[g1]
        @test col isa Dict{Symbol,Vector}
        @test Set(keys(col)) ==
              Set((:topology, :dN, :dS, :id_subject_query, :id_query_subject))
        @test col[:topology] == Vector(pg.topology[:, g1])
        @test col[:dN] == Vector(pg.dN[:, g1])
        @test col[:id_subject_query] == Vector(pg.id_subject_query[:, g1])

        n = size(pg.topology, 1)
        @test_throws ArgumentError pg[n+1]
        @test_throws ArgumentError pg[0]
    end

    @testset "pg[id::AbstractString] - matches integer form" begin
        @test pg["g1"] == pg[g1]
        @test_throws ArgumentError pg["no_such_gene"]
    end

    @testset "pg[indices] - sub-group, matrices sliced correctly" begin
        sub = pg[[g1, g2]]
        @test sub isa PG
        @test sort(collect(keys(sub.intervals))) == ["chr1"]
        @test size(sub.topology) == (2, 2)

        sq, ss = gene_index(sub, "chr1", 100), gene_index(sub, "chr1", 300)
        @test sub.topology[sq, ss] && sub.topology[ss, sq]
        @test sub.dN[min(sq, ss), max(sq, ss)] == pg.dN[min(g1, g2), max(g1, g2)]
        @test sub.id_subject_query[sq, ss] == pg.id_subject_query[g1, g2]
        @test sub.id_query_subject[ss, sq] == pg.id_query_subject[g2, g1]

        # intervals/scaffold_ranges stay concordant in the sub-group too.
        lengths_match = [
            length(rng) == length(sub.intervals[name]) for
            (name, rng) in sub.scaffold_ranges
        ]
        @test all(lengths_match)
    end

    @testset "pg[range] behaves like pg[vector]" begin
        @test pg[1:2].topology == pg[[1, 2]].topology
    end

    @testset "pg[indices] - duplicates and unsorted input collapse" begin
        @test pg[[2, 1, 1]].topology == pg[[1, 2]].topology
    end

    @testset "pg[indices] - bounds and emptiness errors" begin
        n = size(pg.topology, 1)
        @test_throws ArgumentError pg[[0, 1]]
        @test_throws ArgumentError pg[[n+1]]
        @test_throws ArgumentError pg[Int[]]
    end

    @testset "pg[ids::Vector{String}] - matches integer-vector form" begin
        sub_by_id = pg[["g1", "g2"]]
        sub_by_idx = pg[[g1, g2]]
        @test sub_by_id.topology == sub_by_idx.topology
        @test sub_by_id.dN == sub_by_idx.dN

        @test_throws ArgumentError pg[["g1", "ghost"]]
        @test_throws ArgumentError pg[String[]]
    end

    @testset "sub-group with only some relations present" begin
        df2 = DataFrame(q = ["g1"], s = ["g3"], dS = [0.7])
        pg2 = PG(genome, df2)
        sub = pg2[["g1", "g3"]]
        @test sub.dS !== nothing
        @test sub.dN === nothing
        @test sub.id_subject_query === nothing
    end
end

# ============================================================================
# Tests for the *_graph conversion functions
# ============================================================================
@testset "ParalogGroup graph conversions" begin
    PG = BioinfoTools2.Homologs.Paralogs.ParalogGroup
    genes = [
        ("chr1", 100, 200, "+", "g1"),
        ("chr1", 300, 400, "-", "g2"),
        ("chr1", 500, 600, "+", "g3"),
    ]
    genome = genome_from_genes(genes)
    df = DataFrame(
        q = ["g1", "g2"],
        s = ["g2", "g3"],
        dN = [0.1, 0.2],
        dS = [0.3, 0.4],
        id_subject_query = [90.0, 91.0],
        id_query_subject = [80.0, 81.0],
    )
    pg = PG(genome, df)
    g1 = gene_index(pg, "chr1", 100)
    g2 = gene_index(pg, "chr1", 300)
    g3 = gene_index(pg, "chr1", 500)

    @testset "topology_graph - unweighted" begin
        g = topology_graph(pg)
        @test g isa SimpleGraph
        @test nv(g) == 3
        @test ne(g) == 2
        @test has_edge(g, g1, g2) && has_edge(g, g2, g3)
        @test !has_edge(g, g1, g3)
    end

    @testset "dN_graph / dS_graph - symmetric weighted" begin
        dng = dN_graph(pg)
        @test dng isa SimpleWeightedGraph
        @test nv(dng) == 3
        @test get_weight(dng, g1, g2) == 0.1
        @test get_weight(dng, g2, g1) == 0.1   # symmetric both ways
        @test get_weight(dng, g2, g3) == 0.2

        dsg = dS_graph(pg)
        @test dsg isa SimpleWeightedGraph
        @test get_weight(dsg, g1, g2) == 0.3
        @test get_weight(dsg, g3, g2) == 0.4
    end

    @testset "id_*_graph - directed weighted, not mirrored" begin
        sqg = id_subject_query_graph(pg)
        @test sqg isa SimpleWeightedDiGraph
        @test get_weight(sqg, g1, g2) == 90.0
        @test get_weight(sqg, g2, g1) == 0.0   # directed: opposite entry unset

        qsg = id_query_subject_graph(pg)
        @test qsg isa SimpleWeightedDiGraph
        @test get_weight(qsg, g2, g1) == 80.0
        @test get_weight(qsg, g1, g2) == 0.0
    end

    @testset "absent relations propagate as nothing" begin
        pg2 = PG(genome, DataFrame(q = ["g1"], s = ["g2"]))
        @test topology_graph(pg2) isa SimpleGraph
        @test dN_graph(pg2) === nothing
        @test dS_graph(pg2) === nothing
        @test id_subject_query_graph(pg2) === nothing
        @test id_query_subject_graph(pg2) === nothing
    end

    @testset "graphs stay consistent through sub-group indexing" begin
        sub = pg[["g1", "g2"]]
        sg1, sg2 = gene_index(sub, "chr1", 100), gene_index(sub, "chr1", 300)
        @test nv(topology_graph(sub)) == 2
        @test get_weight(dN_graph(sub), sg1, sg2) == 0.1
        @test get_weight(id_subject_query_graph(sub), sg1, sg2) == 90.0
    end
end

# ============================================================================
# Tests for rbh(::ParalogGroup) - per-component reciprocal best hit
# ============================================================================
@testset "rbh(::ParalogGroup)" begin
    PG = BioinfoTools2.Homologs.Paralogs.ParalogGroup
    genes = [
        ("chr1", 100, 110, "+", "g1"),
        ("chr1", 200, 210, "+", "g2"),
        ("chr1", 300, 310, "+", "g3"),
        ("chr1", 400, 410, "+", "g4"),
        ("chr2", 10, 20, "+", "g5"),
        ("chr2", 30, 40, "+", "g6"),
    ]
    genome = genome_from_genes(genes)

    @testset "ranks by dS, then dN, then %ID, per component" begin
        # Component {g1,g2,g3}: g2-g3 has the lowest dS -> wins outright.
        # Component {g4,g5,g6} is a chain g4-g5, g5-g6 with equal dS; dN then
        # breaks the tie in favour of g5-g6.
        df = DataFrame(
            q = ["g1", "g1", "g2", "g4", "g5"],
            s = ["g2", "g3", "g3", "g5", "g6"],
            dS = [0.5, 0.5, 0.1, 0.3, 0.3],
            dN = [0.2, 0.2, 0.2, 0.9, 0.4],
        )
        pg = PG(genome, df)
        result = rbh(pg)

        @test nrow(result) == 2
        @test Set(names(result)) == Set(["query", "subject", "dN", "dS"])

        r1 = only(
            filter(r -> Set([r.query, r.subject]) == Set(["g2", "g3"]), eachrow(result)),
        )
        @test r1.dS == 0.1

        r2 = only(
            filter(r -> Set([r.query, r.subject]) == Set(["g5", "g6"]), eachrow(result)),
        )
        @test r2.dS == 0.3 && r2.dN == 0.4
    end

    @testset "%ID breaks ties after dS and dN, scoring is configurable" begin
        # g1-g2 and g1-g3 tie on both dS and dN; %ID must decide.
        df = DataFrame(
            q = ["g1", "g1"],
            s = ["g2", "g3"],
            dS = [0.2, 0.2],
            dN = [0.1, 0.1],
            id_subject_query = [70.0, 95.0],
            id_query_subject = [72.0, 60.0],
        )
        pg = PG(genome, df)

        # mean(70,72)=71 vs mean(95,60)=77.5 -> g1-g3 wins on higher mean %ID.
        mean_result = rbh(pg; scoring = "mean")
        @test nrow(mean_result) == 1
        @test Set([mean_result.query[1], mean_result.subject[1]]) == Set(["g1", "g3"])

        # min(70,72)=70 vs min(95,60)=60 -> g1-g2 wins on higher minimum %ID.
        min_result = rbh(pg; scoring = "min")
        @test Set([min_result.query[1], min_result.subject[1]]) == Set(["g1", "g2"])

        # max(70,72)=72 vs max(95,60)=95 -> g1-g3 wins on higher maximum %ID.
        max_result = rbh(pg; scoring = "max")
        @test Set([max_result.query[1], max_result.subject[1]]) == Set(["g1", "g3"])

        @test_throws ArgumentError rbh(pg; scoring = "bogus")
    end

    @testset "genuine ties warn and resolve deterministically" begin
        # g1-g2 and g1-g3 are identical on every available metric.
        df = DataFrame(q = ["g1", "g1"], s = ["g2", "g3"], dS = [0.4, 0.4], dN = [0.1, 0.1])
        pg = PG(genome, df)

        local result
        result = @test_logs (:warn, r"tied best-hit") rbh(pg)
        @test nrow(result) == 1
        # The lowest-indexed pair (g1-g2) deterministically wins the tie.
        @test Set([result.query[1], result.subject[1]]) == Set(["g1", "g2"])

        # No tie -> no warning at all.
        clean = PG(genome, DataFrame(q = ["g1"], s = ["g2"], dS = [0.4]))
        @test_logs rbh(clean)
    end

    @testset "singleton (edge-less) components are skipped" begin
        # g4/g5/g6 never appear in the pair table, so they never enter `pg`;
        # only the connected g1-g2 pair contributes a row.
        pg = PG(genome, DataFrame(q = ["g1"], s = ["g2"], dS = [0.1]))
        result = rbh(pg)
        @test nrow(result) == 1
        @test Set([result.query[1], result.subject[1]]) == Set(["g1", "g2"])
    end

    @testset "no rankable relation errors" begin
        pg = PG(genome, DataFrame(q = ["g1"], s = ["g2"]))
        @test_throws ArgumentError rbh(pg)
    end

    @testset "output round-trips through the ParalogGroup constructor" begin
        df = DataFrame(q = ["g1", "g1"], s = ["g2", "g3"], dS = [0.5, 0.1], dN = [0.2, 0.2])
        pg = PG(genome, df)
        result = rbh(pg)
        pg2 = PG(genome, result)
        @test size(pg2.topology, 1) == 2
    end
end
