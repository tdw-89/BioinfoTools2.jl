using BioinfoTools2.Plotting
using Test

# CairoMakie is a package dependency, not a test one; reach it through Plotting.
const Makie = BioinfoTools2.Plotting.CairoMakie

@testset "Plotting" begin
    @testset "metagene_xticks and metagene_boundaries" begin
        positions, labels = metagene_xticks(500, 100)
        @test positions == [1.0, 500.5, 600.5, 1100.0]
        @test labels == ["-500 bp", "TSS", "TES", "+500 bp"]
        @test metagene_boundaries(500, 100) == positions[2:3]
    end

    @testset "stretch_body_columns - flanks kept, body repeated" begin
        # flank 1, body 2 bins stretched to 4 columns.
        matrix = [10 1 2 20; 30 3 4 40]
        stretched = stretch_body_columns(matrix, 1, 2, 4)
        @test stretched == [10 1 1 2 2 20; 30 3 3 4 4 40]
        @test stretch_body_columns(matrix, 1, 2, 2) === matrix
        @test_throws DimensionMismatch stretch_body_columns(matrix, 1, 3, 4)
    end

    @testset "drawing functions place one axis each" begin
        figure = Makie.Figure()
        lines_axis = metagene_lines!(figure[1, 1], rand(7); flank = 2, title = "lines")
        heat_axis, heat =
            metagene_heatmap!(figure[2, 1], [rand(2, 7); fill(NaN, 1, 7)]; flank = 2)
        Makie.Colorbar(figure[2, 2], heat)
        violin_axis = violin_box!(figure[3, 1], [1, 1, 2, 2], [0.1, 0.2, 0.3, 0.4])
        empty_axis = violin_box!(figure[4, 1], Int[], Float64[])

        axes = (lines_axis, heat_axis, violin_axis, empty_axis)
        @test all(axis -> axis isa Makie.Axis, axes)
        @test heat_axis.xticks[] == metagene_xticks(2, 3)
        @test heat_axis.yticks[] == 1:3
        @test violin_axis.xticks[] == 1:2
        # Without boundaries an axis holds only its data plot.
        bare_axis, _ =
            metagene_heatmap!(figure[5, 1], rand(2, 7); flank = 2, mark_boundaries = false)
        bare_lines =
            metagene_lines!(figure[6, 1], rand(7); flank = 2, mark_boundaries = false)
        @test length(bare_axis.scene.plots) == 1
        @test length(bare_lines.scene.plots) == 1
        @test length(heat_axis.scene.plots) == 2

        # Headless render: the whole figure draws without error.
        @test Makie.colorbuffer(figure) isa AbstractMatrix
    end
end
