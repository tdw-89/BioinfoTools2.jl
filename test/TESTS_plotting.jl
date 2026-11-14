using BioinfoTools2.Plotting
using BioinfoTools2.Orthologs
using BioinfoTools2.Reference
using Test

# CairoMakie is a package dependency, not a test one; reach it through Plotting.
const Makie = BioinfoTools2.Plotting.CairoMakie

"""
Given a directory, write a three-species OrthoFinder run (`out` genome-less) and
NCBI-style GFF3s for `ape` and `big`. Returns `(directory, gffs)`.
"""
function family_fixture(directory = mktempdir())
    mkpath(joinpath(directory, "Species_Tree"))
    mkpath(joinpath(directory, "Resolved_Gene_Trees"))
    write(
        joinpath(directory, "Species_Tree", "SpeciesTree_rooted_node_labels.txt"),
        "(out:0.5,(ape:0.1,big:0.1)N1:0.3)N0;",
    )
    write(
        joinpath(directory, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt"),
        "OG0000001: (out_O1.1:0.5,((ape_A1.1:0.1,big_B1.1:0.1)n2:0.1,(ape_A2.1:0.1,big_B2.1:0.1)n3:0.1)n1:0.2)n0;\n",
    )
    gff_lines(prefix, start) = [
        "chr1\tRefSeq\tgene\t$start\t$(start+99)\t.\t+\t.\tID=gene-$prefix;gene_biotype=protein_coding",
        "chr1\tRefSeq\tmRNA\t$start\t$(start+99)\t.\t+\t.\tID=rna-$prefix;Parent=gene-$prefix",
        "chr1\tRefSeq\tCDS\t$start\t$(start+99)\t.\t+\t0\tID=cds-$prefix;Parent=rna-$prefix;protein_id=$prefix.1",
    ]
    gffs = Dict(name => joinpath(directory, "$name.gff3") for name in ("ape", "big"))
    write(
        gffs["ape"],
        join(vcat(gff_lines("A1", 1000), gff_lines("A2", 2000)), '\n') * '\n',
    )
    write(
        gffs["big"],
        join(vcat(gff_lines("B1", 1000), gff_lines("B2", 2000)), '\n') * '\n',
    )
    return directory, gffs
end

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

    @testset "gene family phylogram" begin
        # Root 1 over leaf 2 and node 3, which holds leaves 4 and 5.
        y = phylogram_layout(Int32[5, 2, 5, 4, 5])
        @test y == [1.75, 1.0, 2.5, 2.0, 3.0]
        solid, dashed = phylogram_segments(
            Int32[0, 1, 1, 3, 3],
            Float32[NaN, 1, 3, 1, 1],
            [0.0, 1.0, 3.0, 4.0, 4.0],
            y;
            saturation = 2,
        )
        @test dashed == Makie.Point2f[(0, 2.5), (3, 2.5)] && length(solid) == 10

        directory, gffs = family_fixture()
        clade = load_clade(directory, gffs; orthogroups = nothing)
        tree = divergence_tree(
            clade,
            "OG0000001";
            dS = "((A1.1:0.3,B1.1:0.1):0.2,(A2.1:0.2,B2.1:2.5):0.6,O1.1:1.4);",
        )
        tips = leaves(tree)
        tracks = [
            "5mC" => Dict(
                "ape" => Dict("gene-A1" => 0.8, "gene-A2" => 0.1),
                "big" => Dict("gene-B1" => 0.5),
            ),
        ]
        @test isequal(track_values(tree, tips, tracks), [NaN 0.8 0.5 0.1 NaN])
        @test default_leaf_label.(Ref(tree), tips) == ["O1.1", "A1", "B1", "A2", "B2"]
        palette = species_palette(clade)
        @test length(palette) == 3 && palette["out"] != palette["ape"]

        figure = Makie.Figure(size = (900, 400))
        drawn = gene_family_plot!(figure[1, 1], tree; tracks, layout = :phylogram)
        @test only(drawn.colorbars) isa Makie.Colorbar && drawn.legend isa Makie.Legend
        tree_plots = drawn.tree_axis.scene.plots
        leaf_plot = last(filter(plot -> plot isa Makie.Scatter, tree_plots))
        tooltip = leaf_plot.inspector_label[](leaf_plot, 2, Makie.Point3f(0))
        @test tooltip ==
              "gene-A1\nA1.1 · ape\n5mC: 0.8\nbranch dS 0.3, dN n/a; from root 0.9"
        # B2.1's 2.5 branch is past the default saturation of 2, so it is dashed.
        @test count(plot -> plot isa Makie.LineSegments, tree_plots) == 2
        heat = only(drawn.track_axis.scene.plots)
        @test startswith(
            heat.inspector_label[](heat, (1, 2), Makie.Point3f(0)),
            "5mC: 0.8\nA1",
        )
        @test Makie.colorbuffer(figure) isa AbstractMatrix

        # Without tracks, the legend takes the colorbar's column; pruning draws a stem.
        bare = Makie.Figure()
        stemmed = gene_family_plot!(
            bare[1, 1],
            prune(tree, "ape");
            max_labels = 1,
            layout = :phylogram,
        )
        @test isempty(stemmed.colorbars) && Makie.colorbuffer(bare) isa AbstractMatrix

        # A GeneTrack keeps its own scale, so two scales get two colorbars; the
        # genome-less `out` is looked up by protein, and notes reach the tooltip.
        dosage = GeneTrack(
            "pTriplo",
            Dict("out" => Dict("O1.1" => 0.4), "ape" => Dict("gene-A1" => 0.9));
            colormap = :magma,
            label = "dosage",
            notes = Dict("ape" => Dict("gene-A1" => "own score")),
        )
        @test isequal(track_values(tree, tips, [dosage]), [0.4 0.9 NaN NaN NaN])
        scored = Makie.Figure(size = (900, 400))
        both = gene_family_plot!(
            scored[1, 1],
            tree;
            tracks = [first(tracks), dosage],
            layout = :phylogram,
        )
        @test length(both.colorbars) == 2 &&
              count(plot -> plot isa Makie.Heatmap, both.track_axis.scene.plots) == 2
        scored_leaves =
            last(filter(plot -> plot isa Makie.Scatter, both.tree_axis.scene.plots))
        @test occursin(
            "5mC: 0.8\npTriplo: 0.9 (own score)\n",
            scored_leaves.inspector_label[](scored_leaves, 2, Makie.Point3f(0)),
        )
        @test Makie.colorbuffer(scored) isa AbstractMatrix

        # Opacity fades each cell: A1 at 0.25, others opaque; a missing value stays grey.
        faded = GeneTrack(
            "pTriplo",
            Dict("ape" => Dict("gene-A1" => 0.9, "gene-A2" => 0.1));
            opacity = Dict("ape" => Dict("gene-A1" => 0.25)),
        )
        colors = BioinfoTools2.Plotting._track_colors(
            faded,
            tree,
            tips,
            track_values(tree, tips, [faded])[1, :],
        )
        @test [color.alpha for color in colors] ≈ [1, 0.25, 1, 1, 1] &&
              colors[1] == BioinfoTools2.Plotting.MISSING_COLOR &&
              colors[2].r != colors[4].r

        # A missing value with a coloured reason takes that colour, one without stays
        # grey; the reasons reach the heatmap tooltip and a legend group.
        flagged = GeneTrack(
            "pHaplo",
            Dict("ape" => Dict("gene-A1" => 0.9));
            why_missing = Dict("big" => Dict("gene-B1" => "no dS", "gene-B2" => "other")),
            missing_colors = ["no dS" => :black],
            missing_label = "no ortholog",
        )
        flagged_colors = BioinfoTools2.Plotting._track_colors(
            flagged,
            tree,
            tips,
            track_values(tree, tips, [flagged])[1, :],
        )
        @test flagged_colors[3] == Makie.RGBAf(0, 0, 0, 1) &&
              flagged_colors[5] == flagged_colors[1] == BioinfoTools2.Plotting.MISSING_COLOR
        plain = GeneTrack("5mC", last(first(tracks)))
        @test [
            label for (label, _) in BioinfoTools2.Plotting._missing_entries([plain, flagged, flagged])
        ] == ["no ortholog", "no dS"]
        marked = Makie.Figure(size = (900, 400))
        drawn_flags = gene_family_plot!(marked[1, 1], tree; tracks = [flagged])
        flag_heat = only(drawn_flags.track_axis.scene.plots)
        @test startswith(
            flag_heat.inspector_label[](flag_heat, (1, 3), Makie.Point3f(0)),
            "pHaplo: n/a (no dS)\nB1",
        )
        @test length(drawn_flags.legend.entrygroups[]) == 3 &&
              Makie.colorbuffer(marked) isa AbstractMatrix

        # The chronogram (the default) puts every gene at age 0 and nodes at their
        # mean-path ages, adds the branch-rate colorbar, and keeps the stub of a stem.
        aged = Makie.Figure(size = (900, 400))
        chrono = gene_family_plot!(aged[1, 1], tree; tracks)
        x, y = chrono.positions
        @test all(iszero, x[tips]) && x ≈ node_ages(tree) && chrono.tree_axis.xreversed[]
        @test length(chrono.colorbars) == 2 &&
              last(chrono.colorbars).label[] == "branch rate (relative to the ages)"
        first_tooltips = [
            plot.inspector_label[](plot, 1, Makie.Point3f(0)) for
            plot in chrono.tree_axis.scene.plots if plot isa Makie.Scatter
        ]
        @test any(
            contains(
                "speciation at N0 (OrthoFinder node)\n5 genes below; age (mean dS to them) 1.46",
            ),
            first_tooltips,
        )
        pruned = Makie.Figure()
        @test gene_family_plot!(pruned[1, 1], prune(tree, "ape")).positions[1] ≈
              node_ages(prune(tree, "ape")) && Makie.colorbuffer(aged) isa AbstractMatrix
        @test_throws ArgumentError gene_family_plot!(
            Makie.Figure()[1, 1],
            tree;
            layout = :fan,
        )
        # log2 of length over span, clamped to a quarter and four times.
        @test (
            relative_rate(0.2, 0.1),
            relative_rate(0.1, 0.4),
            relative_rate(1.0, 0.01),
            relative_rate(0.0, 0.0),
            relative_rate(0.3, 0.0),
        ) == (1.0, -2.0, 2.0, 0.0, 2.0)
        @test isnan(relative_rate(NaN, 1.0))
    end
end
