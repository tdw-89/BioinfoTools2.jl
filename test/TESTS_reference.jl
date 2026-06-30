using BioinfoTools2.Reference
using Test

const REF_DATA_DIR = joinpath(@__DIR__, "data")
const GFF_SINGLE  = joinpath(REF_DATA_DIR, "NC_003280.10.gff.gz")  # 4 893 genes, 1 scaffold
const GFF_MULTI   = joinpath(REF_DATA_DIR, "genomic.gff.gz")        # 44 795 genes, 7 scaffolds

@testset "Reference" begin

    # -------------------------------------------------------------------------
    @testset "Species constructor" begin
        s = Species("Caenorhabditis elegans")
        @test s.name      == "Caenorhabditis elegans"
        @test s.taxon_id  == ""
        @test isempty(s.genome.scaffolds)
        @test isempty(s.genome.vocab)
        @test isempty(s.genome.vocab_lookup)
        @test isempty(s.genome.meta_offsets)
        @test isempty(s.genome.meta_blob)

        s2 = Species("Homo sapiens"; taxon_id = "9606")
        @test s2.name     == "Homo sapiens"
        @test s2.taxon_id == "9606"
        @test isempty(s2.genome.scaffolds)
    end

    # -------------------------------------------------------------------------
    @testset "add_features! – single scaffold (NC_003280.10.gff.gz)" begin
        s = Species("C. elegans")
        add_features!(GFF_SINGLE, s.genome)

        # Exactly one scaffold
        @test length(s.genome.scaffolds) == 1
        @test haskey(s.genome.scaffolds, "NC_003280.10")

        # 4 893 gene records parsed:
        #   meta_offsets has one entry per gene + 1 final sentinel
        #   meta_blob has 3 UInt16 tokens (id, source, biotype) × 2 bytes each per gene
        n_genes = 4_893
        @test length(s.genome.meta_offsets) == n_genes + 1
        @test length(s.genome.meta_blob)    == n_genes * 6

        # Vocab is populated (sources and biotypes are interned)
        @test !isempty(s.genome.vocab)

        # The sentinel offset equals (blob length + 1), i.e. points just past the end
        @test s.genome.meta_offsets[end] == length(s.genome.meta_blob) + 1
    end

    # -------------------------------------------------------------------------
    @testset "get_metadata – first gene in NC_003280.10.gff.gz" begin
        # First gene record: gene-CELE_2L52.2, source=RefSeq, biotype=ncRNA
        s = Species("C. elegans")
        add_features!(GFF_SINGLE, s.genome)

        meta = get_metadata(s.genome, UInt32(1))
        @test length(meta) == 3
        @test meta[1] == "gene-CELE_2L52.2"
        @test meta[2] == "RefSeq"
        @test meta[3] == "ncRNA"
    end

    # -------------------------------------------------------------------------
    @testset "get_metadata – repeated calls are consistent" begin
        s = Species("C. elegans")
        add_features!(GFF_SINGLE, s.genome)

        # Calling get_metadata twice for the same index must return equal results
        @test get_metadata(s.genome, UInt32(1)) == get_metadata(s.genome, UInt32(1))
        @test get_metadata(s.genome, UInt32(2)) == get_metadata(s.genome, UInt32(2))

        # Adjacent genes must have individually valid 3-element metadata vectors
        for idx in UInt32(1):UInt32(10)
            m = get_metadata(s.genome, idx)
            @test length(m) == 3
            @test all(!isempty, m)
        end
    end

    # -------------------------------------------------------------------------
    @testset "add_features! – multiple scaffolds (genomic.gff.gz)" begin
        s = Species("C. elegans"; taxon_id = "6239")
        add_features!(GFF_MULTI, s.genome)

        expected_scaffolds = [
            "NC_001328.1", "NC_003279.8", "NC_003280.10",
            "NC_003281.10", "NC_003282.8", "NC_003283.11", "NC_003284.9"
        ]

        @test length(s.genome.scaffolds) == length(expected_scaffolds)
        for id in expected_scaffolds
            @test haskey(s.genome.scaffolds, id)
        end

        # Total gene count across all scaffolds
        n_genes = 44_795
        @test length(s.genome.meta_offsets) == n_genes + 1
        @test length(s.genome.meta_blob)    == n_genes * 6

        # Per-scaffold gene counts (from awk '$3=="gene"' | sort | uniq -c)
        expected_counts = Dict(
            "NC_001328.1"  =>    36,
            "NC_003279.8"  => 4_015,
            "NC_003280.10" => 4_893,
            "NC_003281.10" => 3_681,
            "NC_003282.8"  => 19_323,
            "NC_003283.11" => 7_004,
            "NC_003284.9"  => 5_843,
        )
        total = sum(values(expected_counts))
        @test total == n_genes
    end

    # -------------------------------------------------------------------------
    @testset "add_features! – idempotent genome state after load" begin
        # Loading the same file twice into separate Species must produce
        # identical meta_blob lengths and vocab sizes.
        s1 = Species("C. elegans")
        s2 = Species("C. elegans")
        add_features!(GFF_SINGLE, s1.genome)
        add_features!(GFF_SINGLE, s2.genome)

        @test length(s1.genome.meta_blob)    == length(s2.genome.meta_blob)
        @test length(s1.genome.meta_offsets) == length(s2.genome.meta_offsets)
        @test length(s1.genome.vocab)        == length(s2.genome.vocab)
        @test sort(s1.genome.vocab)          == sort(s2.genome.vocab)
    end

end
