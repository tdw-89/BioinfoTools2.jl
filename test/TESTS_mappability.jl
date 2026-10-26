using BioinfoTools2.Reference
using BioinfoTools2.Data
using BioinfoTools2.Data.Methylation
using BioinfoTools2.Data.Mappability
using BioinfoTools2.Exploration
using IntervalTrees
using StructArrays
using Test
using Zarr

# A UniMap run on the mappability pipeline's test fixture (three chromosomes,
# seven protein-coding genes, 1 kb flanks so that uniqB and nearC overlap), with
# the loci.tsv and GFF3 it was built from. dupA has a paralog on chr2.
const MP_DIR = joinpath(@__DIR__, "data", "unimap")
const MP_ZARR = joinpath(MP_DIR, "mappability.zarr")
const MP_LOCI = joinpath(MP_DIR, "loci.tsv")
const MP_GFF = joinpath(MP_DIR, "annotation.gff3")

"""Load a GFF3 into a fresh genome."""
function mp_genome(path)
    species = Species("fixture")
    add_features!(path, species.genome)
    return species.genome
end

"""Write a copy of the fixture GFF3 with `edit` applied to each gene line."""
function mp_edited_gff(edit)
    path = tempname() * ".gff3"
    open(path, "w") do io
        for line in eachline(MP_GFF)
            fields = split(line, '\t')
            println(io, length(fields) == 9 && fields[3] == "gene" ? edit(fields) : line)
        end
    end
    return path
end

"""Shift a GFF3 feature line by `offset` bases."""
mp_shifted(fields, offset) = join(
    [
        fields[1:3];
        string(parse(Int, fields[4]) + offset);
        string(parse(Int, fields[5]) + offset);
        fields[6:9]
    ],
    '\t',
)

@testset "Mappability" begin
    genome = mp_genome(MP_GFF)
    store = zopen(MP_ZARR; consolidated = true)
    store_loci = store.attrs["loci"]
    read_lengths, fragment_lengths, cutoffs = (
        Int.(store.attrs[key]) for
        key in ("read_lengths", "fragment_lengths", "mapq_cutoffs")
    )
    cache = convert_unimap(MP_ZARR, mktempdir(), genome; channels = :all, loci = MP_LOCI)
    channels = mappability_channels(cache)

    @testset "conversion is lossless" begin
        @test length(channels) ==
              length(read_lengths) * length(fragment_lengths) * length(cutoffs)
        matches = Bool[]
        for channel in channels
            data = read_mappability(cache; channel...)
            iL, iF, ic = findfirst(==(channel.L), read_lengths),
            findfirst(==(channel.F), fragment_lengths),
            findfirst(==(channel.mapq_cutoff), cutoffs)
            for locus in store_loci
                first_base, last_base = locus["start0"] + 1, locus["end0"]
                rows = (locus["offset"]+1):(locus["offset"]+last_base-first_base+1)
                expected = (
                    store["n_simulated"][rows, iF, iL],
                    store["n_correct"][rows, iF, ic, iL],
                    store["n_ectopic"][rows, iF, ic, iL],
                )
                runs = runs_in_range(data, locus["chrom"], first_base, last_base)
                push!(matches, expand_counts(runs, first_base, last_base) == expected)
            end
        end
        @test all(matches)
    end

    @testset "overlapping loci are stored once" begin
        span(id) = (
            locus = only(filter(l -> l["locus_id"] == id, store_loci));
            (locus["start0"], locus["end0"])
        )
        @test span("uniqB")[2] > span("nearC")[1]   # the fixture really overlaps them
        disjoint = Bool[]
        for channel in channels,
            runs in values(read_mappability(cache; channel...).scaffolds)

            push!(
                disjoint,
                issorted(runs.start) && all(runs.start[2:end] .> runs.stop[1:(end-1)]),
            )
        end
        @test all(disjoint)
    end

    @testset "genome matching" begin
        manifest = mappability_manifest(cache)
        @test manifest["genome_match"]["n_matched"] == length(store_loci)
        loose = convert_unimap(MP_ZARR, mktempdir(), genome; channels = [(50, 200, 2)])
        @test occursin(
            "inside locus span",
            mappability_manifest(loose)["genome_match"]["method"],
        )

        no_chr3 = tempname() * ".gff3"
        write(
            no_chr3,
            join(filter(line -> !startswith(line, "chr3"), readlines(MP_GFF)) .* '\n'),
        )
        @test_throws ArgumentError convert_unimap(
            MP_ZARR,
            mktempdir(),
            mp_genome(no_chr3);
            channels = :all,
        )

        shifted = mp_genome(mp_edited_gff(fields -> mp_shifted(fields, 25_000)))
        @test_throws ArgumentError convert_unimap(
            MP_ZARR,
            mktempdir(),
            shifted;
            channels = :all,
            loci = MP_LOCI,
        )
        one_off = mp_genome(
            mp_edited_gff(
                fields ->
                    occursin("ID=gene-dupA", fields[9]) ? mp_shifted(fields, 1) :
                    join(fields, '\t'),
            ),
        )
        @test_logs (:warn, r"1 of 7 loci") match_mode = :any convert_unimap(
            MP_ZARR,
            mktempdir(),
            one_off;
            loci = MP_LOCI,
            channels = [(50, 200, 2)],
        )
    end

    @testset "channels, adding, overwrite and compression" begin
        @test_throws ArgumentError convert_unimap(MP_ZARR, mktempdir(), genome)   # no channels named
        @test_throws ArgumentError convert_unimap(
            MP_ZARR,
            mktempdir(),
            genome;
            channels = [(50, 999, 2)],
        )
        subset_dir = mktempdir()
        convert_unimap(
            MP_ZARR,
            subset_dir,
            genome;
            channels = [(50, 200, 2), (30, 80, 0)],
            compress = nothing,
        )
        @test mappability_channels(subset_dir) ==
              [(L = 50, F = 200, mapq_cutoff = 2), (L = 30, F = 80, mapq_cutoff = 0)]
        @test_throws ArgumentError read_mappability(
            subset_dir;
            L = 50,
            F = 50,
            mapq_cutoff = 2,
        )

        # A second call adds only the channels the cache lacks, leaving the rest.
        first_file = joinpath(subset_dir, "L50_F200_q2", "chr1.arrow")
        written = mtime(first_file)
        convert_unimap(MP_ZARR, subset_dir, genome; channels = [(50, 200, 2), (50, 50, 2)])
        @test length(mappability_channels(subset_dir)) == 3 && mtime(first_file) == written
        @test mappability_manifest(subset_dir)["compress"] == "none"
        added = read_mappability(subset_dir; L = 50, F = 50, mapq_cutoff = 2)
        @test all(
            scaffold ->
                collect(added[scaffold]) ==
                collect(read_mappability(cache; L = 50, F = 50, mapq_cutoff = 2)[scaffold]),
            keys(added),
        )

        # Uncompressed and memory-mapped, or compressed and converted with every
        # other channel: the same runs either way.
        mapped = read_mappability(subset_dir; L = 50, F = 200, mapq_cutoff = 2)
        full = read_mappability(cache; L = 50, F = 200, mapq_cutoff = 2)
        @test all(
            scaffold -> collect(mapped[scaffold]) == collect(full[scaffold]),
            keys(full),
        )

        # A cache of another store is refused unless replaced; a non-cache never is.
        other = mktempdir()
        cp(subset_dir, other; force = true)
        manifest_path = joinpath(other, "manifest.toml")
        write(
            manifest_path,
            replace(
                read(manifest_path, String),
                r"genome_checksum = \"[0-9a-f]*\"" => "genome_checksum = \"other\"",
            ),
        )
        @test_throws ArgumentError convert_unimap(
            MP_ZARR,
            other,
            genome;
            channels = [(30, 50, 0)],
        )
        convert_unimap(MP_ZARR, other, genome; channels = [(30, 50, 0)], overwrite = true)
        @test mappability_channels(other) == [(L = 30, F = 50, mapq_cutoff = 0)]
        not_a_cache = mktempdir()
        touch(joinpath(not_a_cache, "keep.txt"))
        @test_throws ArgumentError convert_unimap(
            MP_ZARR,
            not_a_cache,
            genome;
            channels = :all,
            overwrite = true,
        )
        @test isfile(joinpath(not_a_cache, "keep.txt"))
    end

    data = read_mappability(cache; L = 50, F = 200, mapq_cutoff = 2)

    @testset "queries" begin
        @test isempty(runs_in_range(data, "chr1", 20, 10))
        @test isnothing(runs_in_range(data, "chrX", 1, 10))
        positions = collect(9_000:14_000)   # dupA's locus, and a little either side
        @test isequal(
            mappability_at(data, "chr1", positions),
            per_base(data, "chr1", 9_000, 14_000),
        )
        @test all(isnan, per_base(data, "chr1", 45_000, 45_100))   # no locus there
        @test all(isnan, mappability_at(data, "chrX", [1, 2]))
        run = data["chr1"][1]
        @test map_back(run) == run.n_correct / run.n_simulated
        @test log2_ratio(run, 0.5) ==
              Float32(log2((run.n_correct + 0.5) / (run.n_ectopic + 0.5)))
        @test_throws ArgumentError metric_function(:nonsense, 0.5)
    end

    frequency = feature_frequency(genome, :gene, data; flank = 500)

    @testset "feature_frequency" begin
        @test sort!(collect(keys(frequency.features))) == sort!(
            [
                "gene-" .* ["dupA", "uniqB", "nearC", "edgeD", "revE", "nblkF", "sameG"];
                "gene-pseudo1";
                "rna-lnc1"
            ],
        )
        # Brute force: each region base's genomic position, read back per base.
        agrees = Bool[]
        for (scaffold, tree) in get_feature(genome, :gene), interval in tree
            gene_id = feature_id(genome, interval)
            negative = parse_strand(interval.value) == get_strand('-')
            region_start, region_end = Int(interval.first) - 500, Int(interval.last) + 500
            genomic = negative ? (region_end:-1:region_start) : (region_start:region_end)
            values = [
                base < 1 ? NaN32 : only(mappability_at(data, scaffold, [base])) for
                base in genomic
            ]
            measured = findall(isfinite, values)
            levels = FeatureLevels(frequency.features[gene_id])
            push!(agrees, levels.bases == measured && levels.levels == values[measured])
        end
        @test all(agrees)
    end

    @testset "runs summarise exactly as per-base levels" begin
        per_base_frequency = MethylationFrequency(
            0,
            nothing,
            Dict(id => FeatureLevels(runs) for (id, runs) in frequency.features),
        )
        profiles_agree = [
            isequal(
                gene_profile(runs; flank = 500, body_bins = 20),
                gene_profile(FeatureLevels(runs); flank = 500, body_bins = 20),
            ) for runs in values(frequency.features)
        ]
        @test all(profiles_agree)
        @test isequal(
            mean_gene_profile(frequency; flank = 500, body_bins = 20),
            mean_gene_profile(per_base_frequency; flank = 500, body_bins = 20),
        )
        @test tss_window(frequency; flank = 500, window = 400) ==
              tss_window(per_base_frequency; flank = 500, window = 400)
        @test body_mean(frequency) == body_mean(per_base_frequency)
    end

    @testset "biology: the paralog maps poorly, unique genes don't" begin
        body = body_mean(frequency)
        @test body["gene-dupA"] < 0.5 &&
              body["gene-uniqB"] > 0.99 &&
              body["gene-revE"] > 0.99
        low = mappability_intervals(data, genome; below = 0.5)
        overlaps(scaffold, first_base, last_base) =
            !isempty(
                collect(
                    IntervalTrees.intersect(
                        low.scaffolds[scaffold],
                        (UInt32(first_base), UInt32(last_base)),
                    ),
                ),
            )
        @test overlaps("chr1", 10_500, 12_500) && !overlaps("chr1", 30_001, 33_000)
        intervals_low = [
            map_back(run) < 0.5 for run in data["chr1"] if any(
                interval -> interval.first <= run.start && run.stop <= interval.last,
                low.scaffolds["chr1"],
            )
        ]
        @test all(intervals_low)
    end

    @testset "mask_calls" begin
        # Calls in dupA (low), uniqB (high) and outside any locus, on two scaffolds.
        positions = UInt32[2_000, 11_000, 11_500, 31_000, 45_000]
        payloads = UInt32[pack_payload(3, 1) for _ in positions]
        calls = MethylationData(
            Dict(
                "chr1" => aggregated_calls(positions, payloads),
                "chr9" => aggregated_calls(UInt32[5], UInt32[pack_payload(1, 1)]),
            ),
        )
        values = mappability_at(data, "chr1", positions)
        for (unmeasured, keep_nan) in ((:drop, false), (:keep, true))
            masked = mask_calls(calls, data; below = 0.5, unmeasured)
            expected =
                positions[[isnan(value) ? keep_nan : value >= 0.5 for value in values]]
            @test collect(masked["chr1"].pos) == expected
            @test length(masked["chr9"]) == (keep_nan ? 1 : 0)
        end
    end

    @testset "binned_mappability" begin
        binned = binned_mappability(data; bin_size = 1000)
        # Bin 11 of chr1 is bases 10_001-11_000, inside dupA.
        runs = runs_in_range(data, "chr1", 10_001, 11_000)
        simulated, correct, _ = expand_counts(runs, 10_001, 11_000)
        measured = simulated .> 0
        @test binned["chr1"][11] ≈ sum(correct[measured]) / sum(simulated[measured])
        @test isnan(binned["chr1"][20])   # 19_001-20_000: between loci
        @test length(binned["chr1"]) == 36   # bins stop at the last run (nearC's locus ends at 36_000)
    end
end
