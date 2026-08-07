using BioinfoTools2.Exploration
using BioinfoTools2.Reference
using BioinfoTools2.Data
using BioinfoTools2.Data.Methylation
using DataFrames
using IntervalTrees
using KernelDensity
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
            q = quantiles(tab)

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
            has_all = all(gid -> haskey(byid, gid), gene_ids)
            @test has_all
            correct =
                [byid[gid] == (Float64(i), cld(i, 2)) for (i, gid) in enumerate(gene_ids)]
            @test all(correct)

            # Each quartile bin holds exactly two features.
            bin_counts = [count(t -> t[3] == bin, q) == 2 for bin = 1:4]
            @test all(bin_counts)
        end

        # ---------------------------------------------------------------------
        @testset "custom quantiles count" begin
            q2 = quantiles(tab; quantiles = 2)
            byid = Dict(rec.id => qi for (rec, _, qi) in q2)
            # 1..4 → bin 1, 5..8 → bin 2.
            correct = [byid[gid] == (i <= 4 ? 1 : 2) for (i, gid) in enumerate(gene_ids)]
            @test all(correct)
        end

        # ---------------------------------------------------------------------
        @testset "custom merge function" begin
            # sum of two identical columns == 2i; monotone in i, so the bins are
            # unchanged but the paired value now reflects the merge.
            qs = quantiles(tab; merge = sum)
            byid = Dict(rec.id => (val, qi) for (rec, val, qi) in qs)
            correct =
                [byid[gid] == (Float64(2i), cld(i, 2)) for (i, gid) in enumerate(gene_ids)]
            @test all(correct)
        end

        # ---------------------------------------------------------------------
        @testset "invalid quantiles throws" begin
            @test_throws ArgumentError quantiles(tab; quantiles = 0)
            @test_throws ArgumentError quantiles(tab; quantiles = -3)
        end

        # ---------------------------------------------------------------------
        @testset "no matched samples → empty result" begin
            df_none = DataFrame(sample = ["NO_SUCH_ID"], v1 = [1.0], v2 = [2.0])
            tab_none = load_table(sp.genome, df_none)
            empty_q = quantiles(tab_none)
            @test empty_q isa Vector{Tuple{FeatureRecord,Float64,Int}}
            @test isempty(empty_q)
        end
    end

    # =========================================================================
    @testset "quantiles(data, ranking) - rank-based" begin
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

        # v1 ties the first four rows; v2 (descending within the tie) then
        # decides their relative order. Rows 5-8 are already strictly ordered.
        df = DataFrame(
            sample = vcat(gene_ids, ["NO_SUCH_ID"]),
            v1 = Float64[1, 1, 1, 1, 5, 6, 7, 8, 99],
            v2 = Float64[4, 3, 2, 1, 5, 6, 7, 8, 99],
        )
        tab = load_table(sp.genome, df)

        @testset "ranks by first variable, breaks ties with the next" begin
            q = quantiles(tab, ["v1", "v2"]; quantiles = 4)
            @test q isa Vector{Tuple{FeatureRecord,Int}}
            @test length(q) == 8   # NO_SUCH_ID is unmatched, skipped

            byid = Dict{String,Int}()
            for (rec, qi) in q
                byid[rec.id] = qi
            end

            # Ascending (v1, v2): gene4(1,1) gene3(1,2) gene2(1,3) gene1(1,4) gene5..gene8
            expected_order =
                [gene_ids[4], gene_ids[3], gene_ids[2], gene_ids[1], gene_ids[5:8]...]
            correct = [
                byid[gid] == cld(rank_pos * 4, 8) for
                (rank_pos, gid) in enumerate(expected_order)
            ]
            @test all(correct)
        end

        @testset "unrecognised ranking variable errors" begin
            @test_throws ArgumentError quantiles(tab, ["no_such_variable"])
        end

        @testset "empty ranking errors" begin
            @test_throws ArgumentError quantiles(tab, String[])
        end

        @testset "invalid quantiles count errors" begin
            @test_throws ArgumentError quantiles(tab, ["v1"]; quantiles = 0)
        end

        @testset "no matched samples → empty result" begin
            df_none = DataFrame(sample = ["NO_SUCH_ID"], v1 = [1.0])
            tab_none = load_table(sp.genome, df_none)
            @test isempty(quantiles(tab_none, ["v1"]))
        end
    end

    # =========================================================================
    @testset "quantiles(::DataFrame, ranking) - rank-based" begin
        # dS ties the first four rows; dN then decides their relative order.
        # Rows 5-6 are already strictly ordered by dS alone.
        df = DataFrame(
            query = ["g1", "g2", "g3", "g4", "g5", "g6"],
            subject = ["h1", "h2", "h3", "h4", "h5", "h6"],
            dS = [0.5, 0.5, 0.5, 0.5, 0.1, 0.2],
            dN = [0.4, 0.3, 0.2, 0.1, 0.9, 0.9],
        )

        @testset "ranks by first column, breaks ties with the next" begin
            result = quantiles(df, ["dS", "dN"]; quantiles = 3)
            @test names(result) == vcat(names(df), ["quantile"])
            @test nrow(result) == 6
            @test result.query == df.query   # row order unchanged

            byrow = Dict{String,Int}()
            for row in eachrow(result)
                byrow[row.query] = row.quantile
            end
            # Ascending (dS, dN): g5(0.1) g6(0.2) g4(0.5,0.1) g3(0.5,0.2) g2(0.5,0.3) g1(0.5,0.4)
            expected_order = ["g5", "g6", "g4", "g3", "g2", "g1"]
            correct = [
                byrow[gid] == cld(rank_pos * 3, 6) for
                (rank_pos, gid) in enumerate(expected_order)
            ]
            @test all(correct)
        end

        @testset "input DataFrame is not mutated" begin
            original_names = copy(names(df))
            quantiles(df, ["dS"])
            @test names(df) == original_names
        end

        @testset "remaining ties fall back to original row order" begin
            tied = DataFrame(
                query = ["a", "b", "c"],
                subject = ["x", "y", "z"],
                dS = [0.3, 0.3, 0.3],
            )
            result = quantiles(tied, ["dS"]; quantiles = 3)
            @test result.quantile == [1, 2, 3]   # original order breaks the 3-way tie
        end

        @testset "unrecognised ranking column errors" begin
            @test_throws ArgumentError quantiles(df, ["no_such_column"])
        end

        @testset "empty ranking errors" begin
            @test_throws ArgumentError quantiles(df, String[])
        end

        @testset "invalid quantiles count errors" begin
            @test_throws ArgumentError quantiles(df, ["dS"]; quantiles = 0)
        end

        @testset "empty DataFrame yields an empty (but shaped) result" begin
            empty_df = DataFrame(query = String[], subject = String[], dS = Float64[])
            result = quantiles(empty_df, ["dS"])
            @test nrow(result) == 0
            @test "quantile" in names(result)
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
                cov = coverage(bed_full, :gene)
                @test cov isa Dict{String,Vector{Float64}}
                @test haskey(cov, "NC_003280.10")
                v = cov["NC_003280.10"]
                # Only :gene features can match, so fractions are either 0.0
                # (other feature types) or 1.0 (fully covered genes).
                @test all(x -> x == 0.0 || x == 1.0, v)
                @test count(==(1.0), v) == n_genes
            end

            @testset "filter_zeros drops uncovered features" begin
                vfz = coverage(bed_full, :gene; filter_zeros = true)["NC_003280.10"]
                @test all(==(1.0), vfz)
                @test length(vfz) == n_genes
            end

            @testset "non-overlapping BedData → all zeros" begin
                v = coverage(bed_empty, :gene)["NC_003280.10"]
                @test !isempty(v)
                @test all(==(0.0), v)
            end

            @testset "scaffold absent from BedData is omitted" begin
                cov = coverage(bed_no_match, :gene)
                @test cov isa Dict{String,Vector{Float64}}
                @test isempty(cov)
            end
        end

        @testset "kde" begin
            @testset "non-empty coverage → KDE per scaffold" begin
                k = Exploration.kde(bed_full, :gene)
                @test k isa Dict{String,Union{Nothing,UnivariateKDE}}
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
                k = Exploration.kde(bed_empty, :gene; filter_zeros = true)
                @test haskey(k, "NC_003280.10")
                @test isnothing(k["NC_003280.10"])
            end

            @testset "scaffold absent from BedData is omitted" begin
                @test isempty(Exploration.kde(bed_no_match, :gene))
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

    # =========================================================================
    @testset "methylation" begin
        sp = Species("C. elegans")
        add_features!(EX_GFF_SINGLE, sp.genome)
        const_scaffold = "NC_003280.10"
        scaffold = sp.genome.scaffolds[const_scaffold]

        # Wrap positions/payloads for `const_scaffold` as a MethylationData.
        meth_data(positions, payloads) = MethylationData(
            Dict(const_scaffold => aggregated_calls(UInt32.(positions), payloads)),
        )

        # A stored level sits within half a quantization step (1/255 as a
        # fraction) of the truth, so a sum over `n` sites is within `n`
        # half-steps.
        meth_tol(n_sites) = n_sites * (1 / 255) / 2 + 1e-12

        strand_fwd = BioinfoTools2.BitCodes.STRAND_FWD
        strand_rev = BioinfoTools2.BitCodes.STRAND_REV

        # A gene far enough into the scaffold that a 500 bp left flank fits.
        gene = first(x for x in get_feature(scaffold, :gene) if Int(x.first) > 1000)
        gene_id = Reference.get_metadata_id(sp.genome, Reference.parse_index(gene.value))
        gene_negative = Reference.parse_strand(gene.value) == get_strand('-')
        gene_first, gene_last = Int(gene.first), Int(gene.last)
        gene_length = gene_last - gene_first + 1

        # -------------------------------------------------------------------
        @testset "coverage" begin
            # Two fully methylated sites inside the gene. The score is the sum
            # of per-base fractions over the *feature length*, so uncovered
            # bases pull it down: 2 / gene_length. Same [0, 1] scale as the
            # BedData method.
            data = meth_data(
                [gene_first, gene_first+1],
                UInt32[pack_payload(10, 0), pack_payload(10, 0)],
            )
            scores = coverage(sp.genome, data, :gene)
            @test haskey(scores, const_scaffold)
            @test maximum(scores[const_scaffold]) ≈ 2 / gene_length
            @test all(score -> 0 <= score <= 1, scores[const_scaffold])

            # Half-methylated sites contribute half as much, and a site with no
            # coverage at all (depth 0) contributes nothing. 0.5 is not exactly
            # representable in the 8-bit level field, hence the tolerance.
            half = meth_data(
                [gene_first, gene_first+1, gene_first+2],
                UInt32[pack_payload(5, 5), pack_payload(5, 5), pack_payload(0, 0)],
            )
            half_scores = coverage(sp.genome, half, :gene)
            @test maximum(half_scores[const_scaffold]) ≈ 1 / gene_length atol =
                meth_tol(2) / gene_length

            # Most genes hold no call at all, so most scores are exactly zero.
            @test count(iszero, scores[const_scaffold]) > 0
            filtered = coverage(sp.genome, data, :gene; filter_zeros = true)
            @test all(!iszero, filtered[const_scaffold])
            @test length(filtered[const_scaffold]) < length(scores[const_scaffold])

            # A scaffold the methylation data never mentions is omitted.
            @test isempty(coverage(sp.genome, MethylationData(), :gene))
        end

        # -------------------------------------------------------------------
        @testset "kde" begin
            # Spread sites over several genes so the fitted vector varies.
            positions = Int[]
            payloads = UInt32[]
            for (offset, interval) in enumerate(get_feature(scaffold, :gene))
                offset > 12 && break
                push!(positions, Int(interval.first))
                push!(payloads, pack_payload(offset, 12 - offset))
            end
            data = meth_data(positions, payloads)

            fitted = Exploration.kde(sp.genome, data, :gene)
            @test fitted[const_scaffold] isa UnivariateKDE

            # Nothing left to fit once the zeros are dropped from an empty set.
            @test Exploration.kde(sp.genome, MethylationData(), :gene) |> isempty
        end

        # -------------------------------------------------------------------
        @testset "feature_frequency" begin
            flank = 500
            # One call in the upstream flank, one just inside the 5' end. Both
            # are placed relative to the gene's own orientation.
            upstream = gene_negative ? gene_last + flank : gene_first - flank
            inside = gene_negative ? gene_last - 1 : gene_first + 1

            data = meth_data(
                [upstream, inside],
                UInt32[pack_payload(10, 0), pack_payload(5, 5)],
            )
            frequency = feature_frequency(sp.genome, :gene, data; flank = flank)

            @test frequency isa MethylationFrequency
            @test frequency.min_depth == Exploration.DEFAULT_MIN_DEPTH
            @test frequency.context == CTX_CPG
            @test haskey(frequency.features, gene_id)

            levels = frequency.features[gene_id]
            @test length(levels.levels) == gene_length + 2 * flank
            # Index 1 is the 5' end for either strand.
            @test levels.levels[1] ≈ 1.0
            @test levels.weights[1] == 10
            @test levels.levels[flank+2] ≈ 0.5 atol = 1 / 255
            @test levels.weights[flank+2] == 10

            @testset "depth and context filters" begin
                shallow = meth_data([inside], UInt32[pack_payload(2, 0)])
                @test nnz(
                    feature_frequency(sp.genome, :gene, shallow; flank).features[gene_id].weights,
                ) == 0
                # ...unless the threshold is lowered to admit it.
                @test nnz(
                    feature_frequency(sp.genome, :gene, shallow; flank, min_depth = 2).features[gene_id].weights,
                ) == 1

                non_cpg = meth_data([inside], UInt32[pack_payload(10, 0, CTX_CHH)])
                @test nnz(
                    feature_frequency(sp.genome, :gene, non_cpg; flank).features[gene_id].weights,
                ) == 0
                # `context = nothing` keeps every context.
                kept =
                    feature_frequency(sp.genome, :gene, non_cpg; flank, context = nothing)
                @test nnz(kept.features[gene_id].weights) == 1
                @test kept.context === nothing
            end

            @testset "calls sharing a base combine by depth" begin
                # Same position, opposite strands: 100% over 30 reads and 0%
                # over 10 gives a depth-weighted 0.75 at depth 40.
                shared = meth_data(
                    [inside, inside],
                    UInt32[
                        pack_payload(30, 0, CTX_CPG, strand_fwd),
                        pack_payload(0, 10, CTX_CPG, strand_rev),
                    ],
                )
                combined = feature_frequency(sp.genome, :gene, shared; flank)
                slot = combined.features[gene_id]
                @test slot.weights[flank+2] == 40
                @test slot.levels[flank+2] ≈ 0.75 atol = 1 / 255
            end

            @testset "region is not clipped at the scaffold start" begin
                # A gene closer to the start than `flank`; the region keeps its
                # full width so index `flank + 1` is still the first base.
                edge =
                    first(x for x in get_feature(scaffold, :gene) if Int(x.first) < flank)
                edge_id =
                    Reference.get_metadata_id(sp.genome, Reference.parse_index(edge.value))
                edge_length = Int(edge.last) - Int(edge.first) + 1
                clipped = feature_frequency(sp.genome, :gene, data; flank)
                @test length(clipped.features[edge_id].levels) == edge_length + 2 * flank
            end

            @testset "scaffold absent from the data contributes nothing" begin
                @test isempty(
                    feature_frequency(sp.genome, :gene, MethylationData()).features,
                )
            end
        end

        # -------------------------------------------------------------------
        @testset "gene_profile" begin
            # flank = 2, body = 6 bases over 3 bins → 2 bases per bin.
            #   region:  1 2 | 3 4 5 6 7 8 | 9 10
            #   profile: 1 2 |   3   4   5 | 6  7
            flank, body_bins = 2, 3
            levels = FeatureLevels(
                sparsevec([1, 3, 4, 10], Float32[1.0, 0.25, 0.75, 0.5], 10),
                sparsevec([1, 3, 4, 10], UInt32[10, 5, 15, 8], 10),
            )
            profile = gene_profile(levels; flank, body_bins)

            @test length(profile.levels) == 2 * flank + body_bins
            @test profile.levels[1] ≈ 1.0
            # Bases 3 and 4 share the first body bin, combined by depth:
            # (0.25*5 + 0.75*15) / 20 = 0.625.
            @test profile.levels[3] ≈ 0.625
            @test profile.weights[3] == 20
            # Base 10 is the second downstream flank base → last slot.
            @test profile.levels[7] ≈ 0.5
            @test profile.weights[7] == 8

            # Unmeasured slots carry zero weight, which is what marks them.
            @test profile.weights[2] == 0
            @test all(iszero, profile.weights[[2, 4, 5, 6]])

            @testset "weight transform reshapes the within-bin pooling" begin
                # Bases 3 and 4 again, now pooled under sqrt weighting:
                # (0.25*√5 + 0.75*√15) / (√5 + √15).
                root = gene_profile(levels; flank, body_bins, weight_transform = sqrt)
                expected = (0.25 * sqrt(5) + 0.75 * sqrt(15)) / (sqrt(5) + sqrt(15))
                @test root.levels[3] ≈ expected
                @test root.weights[3] ≈ sqrt(5) + sqrt(15)
                # A single-call slot is untouched by any transform.
                @test root.levels[1] ≈ 1.0

                # Unweighted: both bases count once, so the bin is a plain mean
                # and the weight is the number of measured bases.
                equal = gene_profile(levels; flank, body_bins, weight_by_depth = false)
                @test equal.levels[3] ≈ 0.5
                @test equal.weights[3] == 2
                @test equal.weights[1] == 1
            end

            @testset "invalid weights are rejected" begin
                # log(1) == 0 is allowed (it just drops the base), but a
                # negative or non-finite weight would corrupt the mean.
                @test_throws ArgumentError gene_profile(
                    levels;
                    flank,
                    body_bins,
                    weight_transform = depth -> -depth,
                )
                @test_throws ArgumentError gene_profile(
                    levels;
                    flank,
                    body_bins,
                    weight_transform = depth -> NaN,
                )
                @test Exploration.weight_of(4, sqrt) == 2.0
                @test Exploration.weight_of(1, log) == 0.0
            end

            @testset "too-short regions return nothing" begin
                short = FeatureLevels(spzeros(Float32, 5), spzeros(UInt32, 5))
                @test gene_profile(short; flank, body_bins) === nothing
                exact = FeatureLevels(spzeros(Float32, 6), spzeros(UInt32, 6))
                @test gene_profile(exact; flank, body_bins) !== nothing
            end
        end

        # -------------------------------------------------------------------
        @testset "mean_gene_profile" begin
            flank, body_bins = 2, 3

            # Two genes measured at opposite ends and nowhere else.
            only_start = FeatureLevels(
                sparsevec([1], Float32[1.0], 10),
                sparsevec([1], UInt32[10], 10),
            )
            only_end = FeatureLevels(
                sparsevec([10], Float32[0.5], 10),
                sparsevec([10], UInt32[10], 10),
            )
            frequency = MethylationFrequency(
                Exploration.DEFAULT_MIN_DEPTH,
                CTX_CPG,
                Dict("start" => only_start, "end" => only_end),
            )
            metagene = mean_gene_profile(frequency; flank, body_bins)

            @test length(metagene) == 2 * flank + body_bins
            # Only the gene measured there counts — averaging in a zero for the
            # other would have given 0.5 and 0.25.
            @test metagene[1] ≈ 1.0
            @test metagene[7] ≈ 0.5
            # Nothing measured anywhere else.
            @test all(isnan, metagene[2:6])

            @testset "depth weighting" begin
                # Same base, same 0%/100% split, very different depths. A base
                # measured at 0% still counts: its weight is nonzero even though
                # its level is a structural zero.
                deep = FeatureLevels(
                    sparsevec([1], Float32[1.0], 10),
                    sparsevec([1], UInt32[30], 10),
                )
                shallow = FeatureLevels(
                    sparsevec([1], Float32[0.0], 10),
                    sparsevec([1], UInt32[10], 10),
                )
                pair = MethylationFrequency(
                    Exploration.DEFAULT_MIN_DEPTH,
                    CTX_CPG,
                    Dict("deep" => deep, "shallow" => shallow),
                )

                weighted = mean_gene_profile(pair; flank, body_bins)
                @test weighted[1] ≈ (1.0 * 30 + 0.0 * 10) / 40

                unweighted =
                    mean_gene_profile(pair; flank, body_bins, weight_by_depth = false)
                @test unweighted[1] ≈ 0.5

                # A transform reshapes the cross-gene weighting the same way it
                # reshapes the within-gene pooling: each gene contributes one
                # base here, so its weight is `transform(depth)`.
                compressed =
                    mean_gene_profile(pair; flank, body_bins, weight_transform = sqrt)
                @test compressed[1] ≈
                      (1.0 * sqrt(30) + 0.0 * sqrt(10)) / (sqrt(30) + sqrt(10))
                # sqrt pulls the deep gene's advantage in, so the mean sits
                # between the linear and the equal-weight answers.
                @test unweighted[1] < compressed[1] < weighted[1]

                logged = mean_gene_profile(pair; flank, body_bins, weight_transform = log)
                @test logged[1] ≈ (1.0 * log(30) + 0.0 * log(10)) / (log(30) + log(10))
                @test logged[1] < compressed[1]
            end

            @testset "exclude skips genes" begin
                without_start =
                    mean_gene_profile(frequency; exclude = Set(["start"]), flank, body_bins)
                @test isnan(without_start[1])
                @test without_start[7] ≈ 0.5
            end

            @testset "no qualifying genes → all NaN" begin
                short = FeatureLevels(spzeros(Float32, 5), spzeros(UInt32, 5))
                empty_frequency = MethylationFrequency(
                    Exploration.DEFAULT_MIN_DEPTH,
                    CTX_CPG,
                    Dict("short" => short),
                )
                @test all(isnan, mean_gene_profile(empty_frequency; flank, body_bins))
            end
        end
    end  # methylation
end
