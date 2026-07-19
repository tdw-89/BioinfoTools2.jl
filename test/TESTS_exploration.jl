using BioinfoTools2.Exploration
using BioinfoTools2.Reference
using BioinfoTools2.Data
using DataFrames
using IntervalTrees
using Test

const EX_DATA_DIR = joinpath(@__DIR__, "data")
const EX_GFF_SINGLE = joinpath(EX_DATA_DIR, "NC_003280.10.gff.gz")

@testset "Exploration" begin

    @testset "get_quantiles" begin
        # Load a real genome and grab 8 real gene IDs to build a TabularData with
        # fully predictable per-row means (row i has both columns == i, so the
        # merged mean is exactly i).
        sp = Species("C. elegans")
        add_features!(EX_GFF_SINGLE, sp.genome)
        scaffold = sp.genome.scaffolds["NC_003280.10"]

        gene_ids = String[]
        for iv in get_feature(scaffold, :gene)
            gid = Reference.get_metadata_id(sp.genome, Reference.parse_index(iv.value))
            gid !== nothing && push!(gene_ids, gid)
            length(gene_ids) >= 8 && break
        end
        @test length(gene_ids) == 8

        # 8 matched genes (means 1..8) plus one row that matches nothing.
        df = DataFrame(
            sample = vcat(gene_ids, ["NO_SUCH_ID"]),
            v1 = Float64[1, 2, 3, 4, 5, 6, 7, 8, 99],
            v2 = Float64[1, 2, 3, 4, 5, 6, 7, 8, 99],
        )
        tab = load_table(sp.genome, df, :gene)

        # ---------------------------------------------------------------------
        @testset "default (quantiles = 4, merge = mean)" begin
            q = get_quantiles(sp.genome, tab)

            @test q isa Vector{Tuple{FeatureRecord,Float64,Int}}
            # The unmatched NO_SUCH_ID sample is skipped.
            @test length(q) == 8
            @test !any(t -> t[1].id == "NO_SUCH_ID", q)

            # Every paired record is a resolved gene, every bin is in 1:4.
            @test all(t -> t[1] isa FeatureRecord, q)
            @test all(t -> t[1].feature_type == :gene, q)
            @test all(t -> 1 <= t[3] <= 4, q)

            # For 8 evenly spaced values (1..8) the quartiles pack 2 per bin:
            #   1,2 → q1 | 3,4 → q2 | 5,6 → q3 | 7,8 → q4  (cld(i, 2)).
            byid = Dict(rec.id => (val, qi) for (rec, val, qi) in q)
            for (i, gid) in enumerate(gene_ids)
                @test haskey(byid, gid)
                val, qi = byid[gid]
                @test val == Float64(i)
                @test qi == cld(i, 2)
            end

            # Each quartile bin holds exactly two features.
            for bin = 1:4
                @test count(t -> t[3] == bin, q) == 2
            end
        end

        # ---------------------------------------------------------------------
        @testset "custom quantiles count" begin
            q2 = get_quantiles(sp.genome, tab; quantiles = 2)
            byid = Dict(rec.id => qi for (rec, _, qi) in q2)
            for (i, gid) in enumerate(gene_ids)
                # 1..4 → bin 1, 5..8 → bin 2.
                @test byid[gid] == (i <= 4 ? 1 : 2)
            end
        end

        # ---------------------------------------------------------------------
        @testset "custom merge function" begin
            # sum of two identical columns == 2i; monotone in i, so the bins are
            # unchanged but the paired value now reflects the merge.
            qs = get_quantiles(sp.genome, tab; merge = sum)
            byid = Dict(rec.id => (val, qi) for (rec, val, qi) in qs)
            for (i, gid) in enumerate(gene_ids)
                val, qi = byid[gid]
                @test val == Float64(2i)
                @test qi == cld(i, 2)
            end
        end

        # ---------------------------------------------------------------------
        @testset "invalid quantiles throws" begin
            @test_throws ArgumentError get_quantiles(sp.genome, tab; quantiles = 0)
            @test_throws ArgumentError get_quantiles(sp.genome, tab; quantiles = -3)
        end

        # ---------------------------------------------------------------------
        @testset "no matched samples → empty result" begin
            df_none = DataFrame(sample = ["NO_SUCH_ID"], v1 = [1.0], v2 = [2.0])
            tab_none = load_table(sp.genome, df_none, :gene)
            empty_q = get_quantiles(sp.genome, tab_none)
            @test empty_q isa Vector{Tuple{FeatureRecord,Float64,Int}}
            @test isempty(empty_q)
        end
    end

    # =========================================================================
    @testset "coverage & kde" begin
        sp = Species("C. elegans")
        add_features!(EX_GFF_SINGLE, sp.genome)
        scaffold = sp.genome.scaffolds["NC_003280.10"]
        n_genes = length(get_feature(scaffold, :gene))

        # A BedData whose single interval blankets the whole scaffold, so every
        # gene is fully covered (fraction 1.0).
        bed_full = BedData(
            Dict(
                "NC_003280.10" => let t = IntervalMeta64()
                    push!(t, IntervalValue(UInt32(1), UInt32(100_000_000), UInt64(0)))
                    t
                end,
            ),
        )

        # A BedData on the right scaffold whose interval sits past every feature,
        # so nothing overlaps and all fractions are 0.0.
        bed_empty = BedData(
            Dict(
                "NC_003280.10" => let t = IntervalMeta64()
                    push!(
                        t,
                        IntervalValue(UInt32(200_000_000), UInt32(200_000_001), UInt64(0)),
                    )
                    t
                end,
            ),
        )

        # A BedData whose only scaffold name has no counterpart in the genome.
        bed_no_match = BedData(Dict("NOMATCH" => IntervalMeta64()))

        @testset "coverage" begin
            @testset "full coverage → gene fractions are 1.0" begin
                cov = coverage(sp.genome, :gene, bed_full)
                @test cov isa Dict{String,Vector{Float64}}
                @test haskey(cov, "NC_003280.10")
                v = cov["NC_003280.10"]
                # Only :gene features can match, so fractions are either 0.0
                # (other feature types) or 1.0 (fully covered genes).
                @test all(x -> x == 0.0 || x == 1.0, v)
                @test count(==(1.0), v) == n_genes
            end

            @testset "filter_zeros drops uncovered features" begin
                vfz =
                    coverage(sp.genome, :gene, bed_full; filter_zeros = true)["NC_003280.10"]
                @test all(==(1.0), vfz)
                @test length(vfz) == n_genes
            end

            @testset "non-overlapping BedData → all zeros" begin
                v = coverage(sp.genome, :gene, bed_empty)["NC_003280.10"]
                @test !isempty(v)
                @test all(==(0.0), v)
            end

            @testset "scaffold absent from BedData is omitted" begin
                cov = coverage(sp.genome, :gene, bed_no_match)
                @test cov isa Dict{String,Vector{Float64}}
                @test isempty(cov)
            end
        end

        @testset "kde" begin
            @testset "non-empty coverage → KDE per scaffold" begin
                k = kde(sp.genome, :gene, bed_full)
                @test k isa Dict{String}
                @test haskey(k, "NC_003280.10")
                fit = k["NC_003280.10"]
                # A fitted UnivariateKDE exposes matching grid/density vectors.
                @test !isnothing(fit)
                @test hasproperty(fit, :x) && hasproperty(fit, :density)
                @test length(fit.x) == length(fit.density)
            end

            @testset "empty coverage vector → nothing" begin
                # All fractions are 0.0 and filter_zeros removes them, leaving an
                # empty vector that cannot be fit.
                k = kde(sp.genome, :gene, bed_empty; filter_zeros = true)
                @test haskey(k, "NC_003280.10")
                @test isnothing(k["NC_003280.10"])
            end

            @testset "scaffold absent from BedData is omitted" begin
                @test isempty(kde(sp.genome, :gene, bed_no_match))
            end
        end
    end
end
