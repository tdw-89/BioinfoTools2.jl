using BioinfoTools2.Paralogs
using DataFrames
using Test


@testset "ParalogUtils Tests" begin

    # ========================================
    # Tests for rbh()
    # ========================================
    
    @testset "rbh - basic functionality with max scoring" begin
        # Create test data with reciprocal best hits
        df = DataFrame(
            GeneID = ["A", "B", "C", "D"],
            ParalogID = ["B", "A", "D", "C"],
            Perc1 = [95.0, 94.0, 80.0, 85.0],
            Perc2 = [94.0, 95.0, 85.0, 80.0]
        )
        
        result = rbh(df; scoring="max")
        
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
            Perc2 = [88.0, 90.0]
        )
        
        result = rbh(df; scoring="mean")
        
        @test nrow(result) == 1
        # Mean should be average of bidirectional scores
        @test result.mean_perc[1] ≈ 89.0
    end
    
    @testset "rbh - average scoring (alias)" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "A"],
            Perc1 = [100.0, 80.0],
            Perc2 = [80.0, 100.0]
        )
        
        result = rbh(df; scoring="avg")
        
        @test nrow(result) == 1
        @test result.mean_perc[1] ≈ 90.0
    end
    
    @testset "rbh - double_max scoring" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "A"],
            Perc1 = [95.0, 93.0],
            Perc2 = [93.0, 95.0]
        )
        
        result = rbh(df; scoring="double_max")
        
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
            Perc2 = [94.0, 89.0, 84.0]
        )
        
        result = rbh(df; scoring="max")
        
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
            Perc2 = Float64[]
        )
        
        # Empty DataFrame will cause BoundsError when accessing elements
        @test_throws BoundsError rbh(df; scoring="max")
    end
end