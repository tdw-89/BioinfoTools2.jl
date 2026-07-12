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
        #   meta_blob has 3 UInt32 tokens (id, source, biotype) × 4 bytes each per gene
        n_features = 78_407
        @test length(s.genome.meta_offsets) == n_features + 1
        @test length(s.genome.meta_blob)    == n_features * 12 

        # Vocab is populated (sources and biotypes are interned)
        @test !isempty(s.genome.vocab)

        # The sentinel offset equals (blob length + 1), i.e. points just past the end
        @test s.genome.meta_offsets[end] == length(s.genome.meta_blob) + 1
    end

    # -------------------------------------------------------------------------
    @testset "get_metadata – first gene in NC_003280.10.gff.gz" begin
        # First gene record (second feature): gene-CELE_2L52.2, source=RefSeq, biotype=ncRNA
        s = Species("C. elegans")
        add_features!(GFF_SINGLE, s.genome)

        meta = get_metadata(s.genome, UInt32(2))
        @test length(meta) == 3
        @test meta[1] == "gene-CELE_2L52.2"
        @test meta[2] == "RefSeq"
        @test meta[3] == "ncRNA"
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
        n_features = 531872
        @test length(s.genome.meta_offsets) == n_features + 1
        @test length(s.genome.meta_blob)    == n_features * 12

        # Per-scaffold gene counts (from awk '$3=="gene"' | sort | uniq -c)
        expected_counts = Dict(
            "NC_001328.1"  => 123,
            "NC_003279.8"  => 78_407,
            "NC_003280.10" => 63_135,
            "NC_003281.10" => 118_912,
            "NC_003282.8"  => 108_998,
            "NC_003283.11" => 89_181,
            "NC_003284.9"  => 73_116,
        )
        total = sum(values(expected_counts))
        @test total == n_features
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

    # -------------------------------------------------------------------------
    @testset "get_metadata – all overloads" begin
        s = Species("C. elegans")
        add_features!(GFF_SINGLE, s.genome)
        scaffold  = s.genome.scaffolds["NC_003280.10"]
        n_features = 78_407

        @testset "get_metadata(genome, UInt32) – valid indices" begin
            meta1 = get_metadata(s.genome, UInt32(1))
            @test meta1 isa Vector{String}
            @test length(meta1) == 3   # id, source, biotype

            # Known record verified in the sibling testset above
            @test get_metadata(s.genome, UInt32(2)) == ["gene-CELE_2L52.2", "RefSeq", "ncRNA"]
        end

        @testset "get_metadata(genome, UInt32) – out-of-range indices" begin
            # length(meta_offsets) == n_features + 1; first invalid index is n_features + 1
            @test get_metadata(s.genome, UInt32(n_features + 1))   == String[]
            @test get_metadata(s.genome, UInt32(n_features + 100)) == String[]
        end

        @testset "get_metadata(genome, IntervalMeta64)" begin
            results = get_metadata(s.genome, scaffold.features)
            @test results isa Vector{Vector{String}}
            @test length(results) == n_features
            @test all(r -> length(r) == 3, results)
        end

        @testset "get_metadata(genome, Scaffold) – delegates to IntervalMeta64" begin
            @test get_metadata(s.genome, scaffold) == get_metadata(s.genome, scaffold.features)
        end

        @testset "get_metadata(genome) – full genome" begin
            all_meta = get_metadata(s.genome)
            @test all_meta isa Dict{String, Vector{Vector{String}}}
            @test length(all_meta) == length(s.genome.scaffolds)
            @test haskey(all_meta, "NC_003280.10")
            @test length(all_meta["NC_003280.10"]) == n_features
        end
    end

    # -------------------------------------------------------------------------
    @testset "get_feature – all overloads" begin
        s = Species("C. elegans")
        add_features!(GFF_SINGLE, s.genome)
        scaffold = s.genome.scaffolds["NC_003280.10"]

        @testset "get_feature(scaffold, Symbol) – known term" begin
            tree = get_feature(scaffold, :gene)
            @test !isnothing(tree)
            @test !isempty(tree)
            # Every returned interval must carry the same (gene) SO term code
            gene_code = BioinfoTools2.Reference.parse_so_term(first(tree).value)
            @test all(iv -> BioinfoTools2.Reference.parse_so_term(iv.value) == gene_code, tree)
        end

        @testset "get_feature(scaffold, Symbol) – unknown term" begin
            @test isnothing(get_feature(scaffold, :not_a_real_so_term))
        end

        @testset "get_feature(scaffold, AbstractString) – delegates to Symbol" begin
            @test length(get_feature(scaffold, "gene")) == length(get_feature(scaffold, :gene))
        end

        @testset "get_feature(genome, Symbol) – known term" begin
            result = get_feature(s.genome, :gene)
            @test result isa Dict
            # All scaffolds are always present in the result
            @test length(result) == length(s.genome.scaffolds)
            @test haskey(result, "NC_003280.10")
            @test !isempty(result["NC_003280.10"])
            # Genome-level and scaffold-level must agree on count
            @test length(result["NC_003280.10"]) == length(get_feature(scaffold, :gene))
        end

        @testset "get_feature(genome, Symbol) – unknown term" begin
            @test isempty(get_feature(s.genome, :not_a_real_so_term))
        end

        @testset "get_feature(genome, AbstractString) – delegates to Symbol" begin
            by_symbol = get_feature(s.genome, :gene)
            by_string = get_feature(s.genome, "gene")
            @test keys(by_symbol) == keys(by_string)
            for scaffold_name in keys(by_symbol)
                @test length(by_symbol[scaffold_name]) == length(by_string[scaffold_name])
            end
        end
    end
end
