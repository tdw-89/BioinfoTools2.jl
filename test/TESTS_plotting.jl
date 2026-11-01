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
        asymmetric_positions, asymmetric_labels = metagene_xticks(1000, 100, 500)
        @test asymmetric_positions == [1.0, 1000.5, 1100.5, 1600.0]
        @test asymmetric_labels == ["-1000 bp", "TSS", "TES", "+500 bp"]
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
        asymmetric_axis, _ =
            metagene_heatmap!(figure[5, 2], rand(2, 8); flank = 3, downstream = 1)
        @test asymmetric_axis.xticks[] == metagene_xticks(3, 4, 1)
        @test violin_axis.xticks[] == 1:2
        # Without boundaries an axis holds only its data plot.
        bare_axis, _ =
            metagene_heatmap!(figure[5, 1], rand(2, 7); flank = 2, mark_boundaries = false)
        bare_lines =
            metagene_lines!(figure[6, 1], rand(7); flank = 2, mark_boundaries = false)
        @test length(bare_axis.scene.plots) == 1
        @test length(bare_lines.scene.plots) == 1
        @test length(heat_axis.scene.plots) == 2

        trend_axis =
            trend_plot!(figure[7, 1], 1:4, [1.0, 3.0, 2.0, 4.0], [2.0, 2.0, 3.0, 3.0])
        @test length(trend_axis.scene.plots) == 2    # the points and the trend
        mark_positions!(trend_axis, [1, 3], ["first", "second"])
        @test length(trend_axis.scene.plots) == 4    # plus the rules and their labels
        @test_throws DimensionMismatch mark_positions!(trend_axis, [1, 2], ["only one"])
        @test_throws ArgumentError mark_positions!(trend_axis, [1], ["x"]; side = :top)
        fit_axis = trend_plot!(
            figure[9, 1],
            1:4,
            [1.0, 3.0, 2.0, 4.0],
            [2.0, 2.0, 3.0, 3.0];
            trend_label = "moving average",
        )
        overlay_trend!(fit_axis, 1:4, [1.5, 2.0, 2.5, 3.0]; label = "second trend")
        @test length(fit_axis.scene.plots) == 3
        @test Makie.axislegend(fit_axis) isa Makie.Legend
        left_axis =
            trend_plot!(figure[8, 1], 1:4, [1.0, 3.0, 2.0, 4.0], [2.0, 2.0, 3.0, 3.0])
        mark_positions!(left_axis, [2.5], ["everything left"]; side = :left)
        @test length(left_axis.scene.plots) == 4

        probability_axis = probability_plot!(
            figure[10, 1],
            1:4,
            [0.0, 0.25, 0.75, 1.0],
            [0.1, 0.3, 0.7, 0.9];
            observed_label = "observed rate",
            fitted_label = "logistic fit",
        )
        # The 0.5 rule, the observed rate and the fit.
        @test length(probability_axis.scene.plots) == 3
        @test Makie.axislegend(probability_axis) isa Makie.Legend

        box_axis = quantile_boxplot!(figure[11, 1], [1, 1, 3, 3], [0.1, 0.2, 0.3, 0.4])
        @test box_axis.xticks[] == 1:3
        @test length(box_axis.scene.plots) == 1
        @test quantile_boxplot!(figure[11, 2], Int[], Float64[]; n_bins = 4).xticks[] == 1:4

        # 3 blocks over ranks 1:10, flanks of 2 columns standing for 500 bp each.
        rank_axis, rank_plot = rank_heatmap!(
            figure[12, 1],
            [rand(2, 6); fill(NaN, 1, 6)],
            [1, 4, 8, 11];
            flank = 500,
            flank_columns = 2,
            highclip = :red,
        )
        Makie.Colorbar(figure[12, 2], rank_plot)
        @test rank_axis.xticks[] ==
              ([0.5, 2.5, 4.5, 6.5], ["-500 bp", "TSS", "TES", "+500 bp"])
        @test length(rank_axis.scene.plots) == 2    # the heatmap and the TSS/TES rules
        @test_throws DimensionMismatch rank_heatmap!(
            figure[13, 1],
            rand(3, 6),
            [1, 11];
            flank = 500,
            flank_columns = 2,
        )
        mark_rows!(rank_axis, [3.5, 10.4], ["low", "near the top"])
        @test length(rank_axis.scene.plots) == 5    # plus the rules and two labels
        @test_throws DimensionMismatch mark_rows!(rank_axis, [1, 2], ["only one"])

        # Headless render: the whole figure draws without error.
        @test Makie.colorbuffer(figure) isa AbstractMatrix
    end
end
