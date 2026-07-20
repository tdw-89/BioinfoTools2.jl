using BioinfoTools2.Exploration
using BioinfoTools2.Reference
using BioinfoTools2.Data
using DataFrames
using IntervalTrees
using SparseArrays
using Test

const EX_DATA_DIR = joinpath(@__DIR__, "data")
const EX_GFF_SINGLE = joinpath(EX_DATA_DIR, "NC_003280.10.gff.gz")

@testset "Exploration" begin

    @testset "quantiles" begin
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
        tab = load_table(sp.genome, df)

        # ---------------------------------------------------------------------
        @testset "default (quantiles = 4, merge = mean)" begin
            q = quantiles(sp.genome, tab)

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
            q2 = quantiles(sp.genome, tab; quantiles = 2)
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
            qs = quantiles(sp.genome, tab; merge = sum)
            byid = Dict(rec.id => (val, qi) for (rec, val, qi) in qs)
            for (i, gid) in enumerate(gene_ids)
                val, qi = byid[gid]
                @test val == Float64(2i)
                @test qi == cld(i, 2)
            end
        end

        # ---------------------------------------------------------------------
        @testset "invalid quantiles throws" begin
            @test_throws ArgumentError quantiles(sp.genome, tab; quantiles = 0)
            @test_throws ArgumentError quantiles(sp.genome, tab; quantiles = -3)
        end

        # ---------------------------------------------------------------------
        @testset "no matched samples → empty result" begin
            df_none = DataFrame(sample = ["NO_SUCH_ID"], v1 = [1.0], v2 = [2.0])
            tab_none = load_table(sp.genome, df_none)
            empty_q = quantiles(sp.genome, tab_none)
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
            sp.genome,
            Dict(
                "NC_003280.10" => let t = IntervalTreeM64()
                    push!(t, IntervalValue(UInt32(1), UInt32(100_000_000), UInt64(0)))
                    t
                end,
            ),
        )

        # A BedData on the right scaffold whose interval sits past every feature,
        # so nothing overlaps and all fractions are 0.0.
        bed_empty = BedData(
            sp.genome,
            Dict(
                "NC_003280.10" => let t = IntervalTreeM64()
                    push!(
                        t,
                        IntervalValue(UInt32(200_000_000), UInt32(200_000_001), UInt64(0)),
                    )
                    t
                end,
            ),
        )

        # A BedData whose only scaffold name has no counterpart in the genome.
        bed_no_match = BedData(sp.genome, Dict("NOMATCH" => IntervalTreeM64()))

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

    # =========================================================================
    @testset "calculate_frequency" begin
        test_genome = Species("test").genome

        # Build a BedData from `scaffold => [(start, end), ...]` pairs.
        make_bed(pairs...) = begin
            scaffolds = Dict{String,IntervalTreeM64}()
            for (name, ivs) in pairs
                tree = IntervalTreeM64()
                for (s, e) in ivs
                    push!(tree, IntervalValue(UInt32(s), UInt32(e), UInt64(0)))
                end
                scaffolds[name] = tree
            end
            BedData(test_genome, scaffolds)
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
            measurements = [make_bed("chr1" => [(1, 2)]) for _ = 1:256]
            freq = calculate_frequency(measurements)
            @test eltype(freq["chr1"]) == UInt16
            @test Vector(freq["chr1"]) == UInt16[256, 256]
        end
    end  # calculate_frequency

    # =========================================================================
    @testset "feature_frequency" begin
        sp = Species("C. elegans")
        add_features!(EX_GFF_SINGLE, sp.genome)
        scaffold = sp.genome.scaffolds["NC_003280.10"]
        genes = get_feature(scaffold, :gene)

        flank = 500
        # A gene comfortably past the scaffold start so its left flank is not clipped.
        iv = first(x for x in genes if Int(x.first) > flank)
        gene_id = Reference.get_metadata_id(sp.genome, Reference.parse_index(iv.value))
        negative = Reference.parse_strand(iv.value) == get_strand('-')

        region_start = Int(iv.first) - flank
        region_end = Int(iv.last) + flank
        region_len = region_end - region_start + 1

        # Place a known raw count at each end of the padded region so the result's
        # orientation and values can be checked exactly.
        counts = spzeros(UInt16, region_end + 10)
        counts[region_start] = 3
        counts[region_end] = 7
        frequency = Dict("NC_003280.10" => counts)

        ff = feature_frequency(sp.genome, :gene, frequency, 4; flank = flank)

        @test ff isa FeatureFrequency
        @test ff.n == 4                       # measurement count stored, not divided out
        @test haskey(ff.features, gene_id)

        v = ff.features[gene_id]
        @test eltype(v) == UInt32
        @test length(v) == region_len

        if negative
            # Reversed so index 1 stays at the 5' end.
            @test v[region_len] == 3
            @test v[1] == 7
        else
            @test v[1] == 3
            @test v[region_len] == 7
        end

        @testset "unresolved scaffolds contribute nothing" begin
            empty_ff = feature_frequency(
                sp.genome,
                :gene,
                Dict{String,SparseVector{UInt16,Int}}(),
                2,
            )
            @test empty_ff.n == 2
            @test all(v -> nnz(v) == 0, values(empty_ff.features))
        end
    end  # feature_frequency

    # =========================================================================
    @testset "gene_profile" begin
        # flank = 2, body_bins = 3 keeps the arithmetic checkable by hand.
        counts = [2.0, 4.0, 0.0, 6.0, 12.0, 0.0, 8.0, 10.0]  # length 8
        profile = gene_profile(counts, 2; flank = 2, body_bins = 3)

        @test length(profile) == 2 * 2 + 3
        # freq = counts ./ 2 = [1,2, 0,3,6,0, 4,5]; flanks kept, body [0,3,6,0]
        # interpolated onto 3 points ([0, 0.5, 1]) → [0, 4.5, 0].
        @test profile ≈ [1.0, 2.0, 0.0, 4.5, 0.0, 4.0, 5.0]

        @testset "counts divided by n_measurements" begin
            @test gene_profile(counts, 4; flank = 2, body_bins = 3) ≈ profile ./ 2
        end

        @testset "too-short vectors return nothing" begin
            # Need at least 2*flank + 2 entries to fit a >=2 bp body.
            @test gene_profile(zeros(2 * 2 + 1), 1; flank = 2, body_bins = 3) === nothing
            @test gene_profile(zeros(2 * 2 + 2), 1; flank = 2, body_bins = 3) !== nothing
        end

        @testset "works on sparse count vectors" begin
            sparse_counts = sparsevec([1, 8], [2.0, 10.0], 8)
            sparse_profile = gene_profile(sparse_counts, 2; flank = 2, body_bins = 3)
            @test length(sparse_profile) == 7
            @test sparse_profile[1] == 1.0     # 2 / 2
            @test sparse_profile[end] == 5.0   # 10 / 2
        end
    end  # gene_profile

    # =========================================================================
    @testset "mean_gene_profile" begin
        flank, body_bins = 2, 3
        gene_a = sparsevec([1, 8], UInt32[4, 8], 8)
        gene_b = sparsevec([1, 8], UInt32[8, 4], 8)
        short = spzeros(UInt32, 3)               # too short → always skipped
        ff = FeatureFrequency(2, Dict("a" => gene_a, "b" => gene_b, "short" => short))

        profile_a = gene_profile(gene_a, ff.n; flank, body_bins)
        profile_b = gene_profile(gene_b, ff.n; flank, body_bins)

        mean_all = mean_gene_profile(ff; flank, body_bins)
        @test length(mean_all) == 2 * flank + body_bins
        @test mean_all ≈ (profile_a .+ profile_b) ./ 2

        @testset "exclude skips genes" begin
            only_b = mean_gene_profile(ff; exclude = Set(["a"]), flank, body_bins)
            @test only_b ≈ profile_b
        end

        @testset "no qualifying genes → all zeros" begin
            empty_ff = FeatureFrequency(2, Dict("short" => short))
            @test mean_gene_profile(empty_ff; flank, body_bins) ==
                  zeros(2 * flank + body_bins)
        end
    end  # mean_gene_profile
end
