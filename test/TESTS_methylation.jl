using BioinfoTools2.BitCodes
using BioinfoTools2.Methylation
using Arrow
using CodecZlib
using GFF3
using Random
using StructArrays
using Test

const MT_DATA_DIR = joinpath(@__DIR__, "data")
const MICRO_BISMARK = joinpath(MT_DATA_DIR, "micro_CpG_OT.txt")

# Reference counts for micro_CpG_OT.txt (chromosome 15, 406 records over 76
# sites), computed independently of the loader by scanning the fixture here.
function reference_counts(path)
    counts = Dict{UInt32,Tuple{Int,Int}}()
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue   # Bismark version header
        fields = split(line, '\t')
        pos = parse(UInt32, fields[4])
        meth, unmeth = get(counts, pos, (0, 0))
        counts[pos] = fields[5] == "Z" ? (meth + 1, unmeth) : (meth, unmeth + 1)
    end
    return counts
end

@testset "Methylation" begin

    # -------------------------------------------------------------------------
    @testset "payload encoding / decoding" begin
        payload = pack_payload(10, 5, CTX_CPG, STRAND_FWD)
        @test get_meth(payload) == 10
        @test get_unmeth(payload) == 5
        @test get_depth(payload) == 15
        @test get_context(payload) == CTX_CPG
        @test context_label(payload) == :CpG
        @test is_forward(payload)
        @test get_strand(payload) == GFF3.GenomicFeatures.STRAND_POS
        @test meth_fraction(payload) ≈ 10 / 15

        # An AggregatedCall is exactly 8 bytes: two UInt32s, no padding.
        @test sizeof(AggregatedCall) == 8
        @test isbitstype(AggregatedCall)
    end

    # -------------------------------------------------------------------------
    @testset "encoding edge cases" begin
        # Zero counts
        zero_payload = pack_payload(0, 0, CTX_CPG, STRAND_FWD)
        @test get_meth(zero_payload) == 0
        @test get_unmeth(zero_payload) == 0
        @test get_depth(zero_payload) == 0
        @test isnan(meth_fraction(zero_payload))

        # Exactly at the 12-bit maximum
        max_payload = pack_payload(4095, 4095, CTX_CHH, STRAND_REV)
        @test get_meth(max_payload) == 4095
        @test get_unmeth(max_payload) == 4095
        @test get_depth(max_payload) == 8190
        @test get_context(max_payload) == CTX_CHH
        @test !is_forward(max_payload)

        # Overflow clamps instead of wrapping into the neighbouring field
        for (meth, unmeth) in ((4096, 0), (0, 4096), (100_000, 999_999), (typemax(Int), 7))
            clamped = pack_payload(meth, unmeth, CTX_CHG, STRAND_BOTH)
            @test get_meth(clamped) == min(meth, MAX_COUNT)
            @test get_unmeth(clamped) == min(unmeth, MAX_COUNT)
            # The fields must not have bled into context/strand.
            @test get_context(clamped) == CTX_CHG
            @test Methylation.get_strand_code(clamped) == STRAND_BOTH
        end

        # Negative counts floor at zero rather than wrapping to 4095
        @test get_meth(pack_payload(-5, 3, CTX_CPG, STRAND_FWD)) == 0
        @test get_unmeth(pack_payload(-5, 3, CTX_CPG, STRAND_FWD)) == 3

        # Every context / strand combination round-trips, and fields are
        # independent of one another.
        for context in (CTX_CPG, CTX_CHG, CTX_CHH, CTX_UNKNOWN),
            strand in (STRAND_FWD, STRAND_REV, STRAND_BOTH, STRAND_NA)

            payload = pack_payload(1234, 2345, context, strand)
            @test get_meth(payload) == 1234
            @test get_unmeth(payload) == 2345
            @test get_context(payload) == context
            @test Methylation.get_strand_code(payload) == strand
            @test get_strand(payload) == decode_strand(strand)
            @test is_forward(payload) == (strand == STRAND_FWD)
            # Bits 28-31 are reserved and must stay clear.
            @test (payload >> 28) == 0
        end
    end

    # -------------------------------------------------------------------------
    @testset "AggregatedCall constructor and accessors" begin
        call = AggregatedCall(1_000_000, 3, 9, CTX_CHG, STRAND_REV)
        @test call.pos == 1_000_000
        @test get_meth(call) == 3
        @test get_unmeth(call) == 9
        @test get_depth(call) == 12
        @test get_context(call) == CTX_CHG
        @test context_label(call) == :CHG
        @test !is_forward(call)
        @test get_strand(call) == GFF3.GenomicFeatures.STRAND_NEG
        @test meth_fraction(call) ≈ 0.25

        # Positions near the top of the UInt32 range are representable
        @test AggregatedCall(typemax(UInt32), 1, 1).pos == typemax(UInt32)
    end

    # -------------------------------------------------------------------------
    @testset "StructArray construction" begin
        calls = alloc_calls(10)
        @test calls isa StructArray{AggregatedCall}
        @test length(calls) == 10
        @test calls.pos isa Vector{UInt32}
        @test calls.payload isa Vector{UInt32}

        @test length(alloc_calls()) == 0

        positions = UInt32[10, 20, 30]
        payloads = [pack_payload(i, i, CTX_CPG, STRAND_FWD) for i = 1:3]
        wrapped = aggregated_calls(positions, payloads)
        @test length(wrapped) == 3
        @test wrapped[2] == AggregatedCall(20, 2, 2, CTX_CPG, STRAND_FWD)
        # Wrapping is a view of the same columns, not a copy.
        @test wrapped.pos === positions

        @test_throws DimensionMismatch aggregated_calls(positions, payloads[1:2])
    end

    # -------------------------------------------------------------------------
    @testset "find_calls_in_range" begin
        positions = UInt32.(collect(100:100:1000))
        payloads = [pack_payload(i, 0, CTX_CPG, STRAND_FWD) for i = 1:10]
        calls = aggregated_calls(positions, payloads)

        # Inclusive on both ends
        hit = find_calls_in_range(calls, 300, 500)
        @test [c.pos for c in hit] == UInt32[300, 400, 500]

        # Bounds that fall between sites
        @test [c.pos for c in find_calls_in_range(calls, 301, 499)] == UInt32[400]

        # Whole array, and beyond either end
        @test length(find_calls_in_range(calls, 0, 100_000)) == 10
        @test isempty(find_calls_in_range(calls, 1, 99))
        @test isempty(find_calls_in_range(calls, 1001, 2000))

        # Single-site and degenerate ranges
        @test length(find_calls_in_range(calls, 700, 700)) == 1
        @test isempty(find_calls_in_range(calls, 700, 699))

        # The result is still a StructArray, but one whose columns are views
        # into the original (no copying), and it decodes correctly.
        @test hit isa StructArray{AggregatedCall}
        @test hit.pos isa SubArray
        @test get_meth(hit[1]) == 3
        # A view can itself be range-queried.
        @test [c.pos for c in find_calls_in_range(hit, 400, 500)] == UInt32[400, 500]

        @test isempty(find_calls_in_range(alloc_calls(0), 1, 10))
    end

    # -------------------------------------------------------------------------
    @testset "find_calls_in_range - millions of positions" begin
        n = 2_000_000
        # Sites every 3bp, so expected hit counts are known exactly.
        positions = UInt32.(1:3:(3n-2))
        payloads = fill(pack_payload(4, 1, CTX_CPG, STRAND_FWD), n)
        calls = aggregated_calls(positions, payloads)

        @test length(find_calls_in_range(calls, 1, 3n)) == n
        @test length(find_calls_in_range(calls, 1, 300)) == 100
        @test length(find_calls_in_range(calls, 4, 300)) == 99

        rng = Random.MersenneTwister(20260727)
        for _ = 1:200
            lo = rand(rng, 1:(3n))
            hi = lo + rand(rng, 0:10_000)
            hits = find_calls_in_range(calls, lo, hi)
            # Sites sit at 3k + 1 for k in 0:n-1, so the expected hit count is
            # just how many k fall in the window.
            first_k = max(0, cld(lo - 1, 3))
            last_k = min(n - 1, fld(hi - 1, 3))
            @test length(hits) == max(0, last_k - first_k + 1)
            @test all(lo <= c.pos <= hi for c in hits)
        end
    end

    # -------------------------------------------------------------------------
    @testset "infer_strand" begin
        @test infer_strand("CpG_OT_sample_bismark_bt2_pe.deduplicated.txt") == STRAND_FWD
        @test infer_strand("CpG_CTOT_sample.txt") == STRAND_FWD
        @test infer_strand("CHH_OB_sample.txt") == STRAND_REV
        @test infer_strand("CHG_CTOB_sample.txt") == STRAND_REV
        @test infer_strand("/some/dir/CpG_OT_sample.txt.gz") == STRAND_FWD
        @test infer_strand("unlabelled_calls.txt") == STRAND_NA
    end

    # -------------------------------------------------------------------------
    @testset "load_bismark - micro_CpG_OT.txt" begin
        data = load_bismark(MICRO_BISMARK)

        @test data isa MethylationData
        @test length(data) == 1
        @test haskey(data, "15")
        @test collect(keys(data)) == ["15"]

        calls = data["15"]
        @test calls isa StructArray{AggregatedCall}
        @test length(calls) == 76
        @test n_sites(data) == 76

        # Positions must come out sorted (the fixture is in read order).
        @test issorted(calls.pos)
        @test allunique(calls.pos)

        # Total depth equals the number of records in the file.
        @test sum(Int(get_depth(c)) for c in calls) == 406

        # Strand comes from the "CpG_OT" file name, context from the call letter.
        @test all(is_forward, calls)
        @test all(c -> get_context(c) == CTX_CPG, calls)

        # Per-site counts match an independent tally of the fixture.
        expected = reference_counts(MICRO_BISMARK)
        @test length(expected) == 76
        for call in calls
            meth, unmeth = expected[call.pos]
            @test get_meth(call) == meth
            @test get_unmeth(call) == unmeth
        end

        # Range queries against real coordinates
        window = find_calls_in_range(data, "15", 31971000, 31971100)
        @test [Int(c.pos) for c in window] == [31971008, 31971022, 31971030, 31971047, 31971069, 31971077]
        @test find_calls_in_range(data, "no_such_scaffold", 1, 100) === nothing
    end

    # -------------------------------------------------------------------------
    @testset "load_bismark - parsing rules" begin
        # All four contexts, both methylation states, from an in-memory file.
        text = """
        Bismark methylation extractor version v0.25.1
        read1\t+\tchr1\t100\tZ
        read2\t-\tchr1\t100\tz
        read3\t+\tchr1\t100\tZ
        read4\t+\tchr1\t200\tX
        read5\t-\tchr1\t200\tx
        read6\t+\tchr2\t50\tH
        read7\t-\tchr2\t50\th
        read8\t+\tchr2\t60\tU
        read9\t-\tchr2\t60\tu
        """
        data = load_bismark(IOBuffer(text); strand = STRAND_REV)

        @test sort(collect(keys(data))) == ["chr1", "chr2"]
        chr1 = data["chr1"]
        @test length(chr1) == 2
        @test chr1[1] == AggregatedCall(100, 2, 1, CTX_CPG, STRAND_REV)
        @test chr1[2] == AggregatedCall(200, 1, 1, CTX_CHG, STRAND_REV)

        chr2 = data["chr2"]
        @test get_context(chr2[1]) == CTX_CHH
        @test get_context(chr2[2]) == CTX_UNKNOWN

        # The same position in two contexts stays two entries, ordered by
        # position, and both are still reachable by a range query.
        mixed = load_bismark(
            IOBuffer("r1\t+\tc\t10\tZ\nr2\t+\tc\t10\tH\nr3\t+\tc\t10\tZ\n");
            strand = STRAND_FWD,
        )
        entries = mixed["c"]
        @test length(entries) == 2
        @test all(e -> e.pos == 10, entries)
        @test sort([get_context(e) for e in entries]) == [CTX_CPG, CTX_CHH]
        @test length(find_calls_in_range(entries, 10, 10)) == 2

        # Header, blank lines, short lines and unknown call letters are skipped.
        messy = load_bismark(
            IOBuffer(
                "Bismark methylation extractor version v0.25.1\n" *
                "\n" *
                "truncated\tline\n" *
                "r1\t+\tc\t5\tQ\n" *          # unrecognised call letter
                "r2\t+\tc\tnotanumber\tZ\n" * # unparseable position
                "r3\t+\tc\t7\tZ\n",
            ),
        )
        @test n_sites(messy) == 1
        @test messy["c"][1] == AggregatedCall(7, 1, 0, CTX_CPG, STRAND_NA)

        # A file with no usable records yields an empty dataset.
        @test n_sites(load_bismark(IOBuffer("Bismark methylation extractor v0.25.1\n"))) ==
              0
    end

    # -------------------------------------------------------------------------
    @testset "load_bismark - gzipped input" begin
        plain = load_bismark(MICRO_BISMARK)

        gz_path = joinpath(mktempdir(), "CpG_OT_micro.txt.gz")
        open(GzipCompressorStream, gz_path, "w") do io
            write(io, read(MICRO_BISMARK))
        end

        gzipped = load_bismark(gz_path)
        @test n_sites(gzipped) == n_sites(plain)
        # The strand is still inferred from the name, ".gz" suffix and all.
        @test all(is_forward, gzipped["15"])
        @test all(gzipped["15"][i] == plain["15"][i] for i in eachindex(plain["15"]))
    end

    # -------------------------------------------------------------------------
    @testset "load_bismark - saturating aggregation" begin
        # More than 4095 reads at one site must clamp, not wrap.
        io = IOBuffer()
        println(io, "Bismark methylation extractor version v0.25.1")
        for i = 1:5000
            println(io, "read$i\t+\tchr1\t42\tZ")
        end
        seekstart(io)
        data = load_bismark(io; strand = STRAND_FWD)

        call = data["chr1"][1]
        @test get_meth(call) == MAX_COUNT
        @test get_unmeth(call) == 0
        @test get_context(call) == CTX_CPG
        @test is_forward(call)
    end

    # -------------------------------------------------------------------------
    @testset "load_bismark - multiple files merged" begin
        # Reuse the fixture as a second, reverse-strand file: same sites, so the
        # merge must keep the two strands as separate entries.
        dir = mktempdir()
        ob_path = joinpath(dir, "CpG_OB_micro.txt")
        cp(MICRO_BISMARK, ob_path)

        merged = load_bismark([MICRO_BISMARK, ob_path])
        @test n_sites(merged) == 152          # 76 sites × 2 strands
        calls = merged["15"]
        @test issorted(calls.pos)
        @test count(is_forward, calls) == 76
        @test count(!is_forward, calls) == 76
        @test sum(Int(get_depth(c)) for c in calls) == 812

        # Merging a file with itself (same strand) sums the counts instead.
        doubled = merge_calls([load_bismark(MICRO_BISMARK), load_bismark(MICRO_BISMARK)])
        single = load_bismark(MICRO_BISMARK)
        @test n_sites(doubled) == 76
        for (a, b) in zip(doubled["15"], single["15"])
            @test a.pos == b.pos
            @test get_meth(a) == 2 * get_meth(b)
            @test get_unmeth(a) == 2 * get_unmeth(b)
        end

        @test n_sites(load_bismark(String[])) == 0
    end

    # -------------------------------------------------------------------------
    @testset "Arrow round-trip" begin
        data = load_bismark(MICRO_BISMARK)
        calls = data["15"]
        dir = mktempdir()

        for compress in (:zstd, :lz4, nothing)
            path = joinpath(dir, "calls_$(something(compress, :none)).arrow")
            @test write_methylation_arrow(path, calls; compress = compress) == path
            @test isfile(path)

            restored = read_methylation_arrow(path)
            @test restored isa StructArray{AggregatedCall}
            @test length(restored) == length(calls)
            @test collect(restored.pos) == collect(calls.pos)
            @test collect(restored.payload) == collect(calls.payload)
            @test all(restored[i] == calls[i] for i in eachindex(calls))

            # Decoded values survive the round-trip, not just the raw bits.
            for (a, b) in zip(restored, calls)
                @test get_meth(a) == get_meth(b)
                @test get_unmeth(a) == get_unmeth(b)
                @test get_context(a) == get_context(b)
                @test Methylation.get_strand_code(a) == Methylation.get_strand_code(b)
            end

            # And binary search still works on the mapped columns.
            @test [c.pos for c in find_calls_in_range(restored, 31971000, 31971100)] == [c.pos for c in find_calls_in_range(calls, 31971000, 31971100)]
        end

        # An empty array is still a valid (readable) file.
        empty_path = joinpath(dir, "empty.arrow")
        write_methylation_arrow(empty_path, alloc_calls(0))
        @test length(read_methylation_arrow(empty_path)) == 0

        # A file that isn't methylation data is rejected rather than mis-read.
        other_path = joinpath(dir, "other.arrow")
        Arrow.write(other_path, (a = [1, 2, 3], b = [4, 5, 6]))
        @test_throws ArgumentError read_methylation_arrow(other_path)
    end

    # -------------------------------------------------------------------------
    @testset "Arrow round-trip - large, zstd, memory-mapped" begin
        # Large enough to exercise multi-batch writing and real compression.
        n = 1_000_000
        rng = Random.MersenneTwister(89)
        positions = UInt32.(cumsum(rand(rng, 1:20, n)))
        payloads = [
            pack_payload(rand(rng, 0:60), rand(rng, 0:60), CTX_CPG, STRAND_FWD) for _ = 1:n
        ]
        calls = aggregated_calls(positions, payloads)

        path = joinpath(mktempdir(), "big.arrow")
        write_methylation_arrow(path, calls; compress = :zstd)
        restored = read_methylation_arrow(path)

        @test length(restored) == n
        @test collect(restored.pos) == positions
        @test collect(restored.payload) == payloads

        # Zstd should comfortably beat the 8-bytes-per-site in-memory footprint.
        @test filesize(path) < 8 * n

        lo, hi = positions[400_000], positions[600_000]
        @test length(find_calls_in_range(restored, lo, hi)) ==
              length(find_calls_in_range(calls, lo, hi))
    end

    # -------------------------------------------------------------------------
    @testset "dataset round-trip (one file per scaffold)" begin
        data = load_bismark(
            IOBuffer(
                "r1\t+\tchr1\t10\tZ\nr2\t+\tchr1\t20\tz\n" *
                "r3\t+\tchr2\t30\tZ\nr4\t+\tscaffold|weird:name\t40\tZ\n",
            );
            strand = STRAND_FWD,
        )

        dir = joinpath(mktempdir(), "methylation")
        @test write_methylation(dir, data) == dir
        @test isfile(joinpath(dir, "scaffolds.tsv"))

        restored = read_methylation(dir)
        @test sort(collect(keys(restored))) == sort(collect(keys(data)))
        @test n_sites(restored) == n_sites(data)
        for chrom in keys(data)
            original, back = data[chrom], restored[chrom]
            @test length(back) == length(original)
            @test all(back[i] == original[i] for i in eachindex(original))
        end

        # Scaffold names that aren't valid file names survive via the manifest.
        @test haskey(restored, "scaffold|weird:name")

        @test_throws ArgumentError read_methylation(mktempdir())
    end
end
