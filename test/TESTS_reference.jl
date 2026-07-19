using BioinfoTools2.Reference
using BioGenerics
using GFF3
using Test

const REF_DATA_DIR = joinpath(@__DIR__, "data")
const GFF_SINGLE = joinpath(REF_DATA_DIR, "NC_003280.10.gff.gz")  # 4 893 genes, 1 scaffold
const GFF_MULTI = joinpath(REF_DATA_DIR, "genomic.gff.gz")        # 44 795 genes, 7 scaffolds

function read_test_gff_record(line::AbstractString)
    gff = tempname() * ".gff3"
    record = GFF3.Record()
    found = false
    try
        open(gff, "w") do io
            write(io, "##gff-version 3\n")
            write(io, line)
            endswith(line, '\n') || write(io, '\n')
        end

        open(gff) do io
            rdr = GFF3.Reader(io)
            try
                while !eof(rdr)
                    read!(rdr, record)
                    if BioGenerics.isfilled(record)
                        found = true
                        break
                    end
                end
            finally
                close(rdr)
            end
        end
    finally
        rm(gff; force = true)
    end

    found || error("No GFF3 record read from test line")
    return record
end

@testset "Reference" begin

    # -------------------------------------------------------------------------
    @testset "Species constructor" begin
        s = Species("Caenorhabditis elegans")
        @test s.name == "Caenorhabditis elegans"
        @test s.taxon_id == ""
        @test isempty(s.genome.scaffolds)
        @test isempty(s.genome.vocab)
        @test isempty(s.genome.vocab_lookup)
        @test isempty(s.genome.meta_offsets)
        @test isempty(s.genome.meta_blob)

        s2 = Species("Homo sapiens"; taxon_id = "9606")
        @test s2.name == "Homo sapiens"
        @test s2.taxon_id == "9606"
        @test isempty(s2.genome.scaffolds)
    end

    # -------------------------------------------------------------------------
    @testset "sanitize_id" begin
        @test Reference.sanitize_id("gene:WBGene00000001") == "WBGene00000001"
        @test Reference.sanitize_id("gene:transcript:TX0001") == "TX0001"
        @test Reference.sanitize_id("mRNA:RNA:CDS:Protein:ABC123") == "ABC123"
        @test Reference.sanitize_id("gene-family:ABC123") == "gene-family:ABC123"
        @test Reference.sanitize_id("ABC123") == "ABC123"
    end

    # -------------------------------------------------------------------------
    @testset "parse_record" begin
        @testset "known SO term with sanitized ID and full metadata" begin
            record = read_test_gff_record(
                "chr1\tRefSeq\tgene\t5\t20\t.\t+\t.\tID=gene:transcript:GENE0001;gene_biotype=protein_coding",
            )

            result = Reference.parse_record(record, UInt32(7))

            @test result isa Reference.ParseResult
            @test result.scaffold_id == "chr1"
            @test result.start_pos == UInt32(5)
            @test result.end_pos == UInt32(20)
            @test Reference.parse_index(result.code) == UInt32(7)
            @test Reference.parse_strand(result.code) ==
                  GFF3.GenomicFeatures.STRAND_POS
            @test Reference.parse_so_term(result.code) == Reference.convert_so_term("gene")
            @test result.id == "GENE0001"
            @test result.source == "RefSeq"
            @test result.biotype == "protein_coding"
        end

        @testset "sanitization disabled preserves prefixed ID" begin
            record = read_test_gff_record(
                "chr1\tRefSeq\tmRNA\t10\t30\t.\t-\t.\tID=transcript:TX0001;gene_biotype=ncRNA",
            )

            result = Reference.parse_record(
                record,
                UInt32(8);
                sanitize_ids = false,
            )

            @test result.id == "transcript:TX0001"
            @test Reference.parse_strand(result.code) == GFF3.GenomicFeatures.STRAND_NEG
        end

        @testset "missing optional fields fall back to NA" begin
            record = read_test_gff_record("chr2\t.\tgene\t1\t9\t.\t.\t.\tName=no_id")

            result = Reference.parse_record(record, UInt32(9))

            @test result.id == "NA"
            @test result.source == "NA"
            @test result.biotype == "NA"
            @test Reference.parse_strand(result.code) == GFF3.GenomicFeatures.STRAND_BOTH
        end

        @testset "multi-valued ID and biotype fall back to NA" begin
            record = read_test_gff_record(
                "chr3\tRefSeq\tgene\t2\t8\t.\t+\t.\tID=gene:ONE,gene:TWO;gene_biotype=type1,type2",
            )

            result = Reference.parse_record(record, UInt32(10))

            @test result.id == "NA"
            @test result.source == "RefSeq"
            @test result.biotype == "NA"
        end

        @testset "unknown SO term is skipped" begin
            record = read_test_gff_record(
                "chr4\tRefSeq\tnot_a_real_so_term\t3\t6\t.\t+\t.\tID=gene:SKIPME;gene_biotype=protein_coding",
            )

            @test isnothing(Reference.parse_record(record, UInt32(11)))
        end
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
        @test length(s.genome.meta_blob) == n_features * 12

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
            "NC_001328.1",
            "NC_003279.8",
            "NC_003280.10",
            "NC_003281.10",
            "NC_003282.8",
            "NC_003283.11",
            "NC_003284.9",
        ]

        @test length(s.genome.scaffolds) == length(expected_scaffolds)
        for id in expected_scaffolds
            @test haskey(s.genome.scaffolds, id)
        end

        # Total gene count across all scaffolds
        n_features = 531872
        @test length(s.genome.meta_offsets) == n_features + 1
        @test length(s.genome.meta_blob) == n_features * 12

        # Per-scaffold gene counts (from awk '$3=="gene"' | sort | uniq -c)
        expected_counts = Dict(
            "NC_001328.1" => 123,
            "NC_003279.8" => 78_407,
            "NC_003280.10" => 63_135,
            "NC_003281.10" => 118_912,
            "NC_003282.8" => 108_998,
            "NC_003283.11" => 89_181,
            "NC_003284.9" => 73_116,
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

        @test length(s1.genome.meta_blob) == length(s2.genome.meta_blob)
        @test length(s1.genome.meta_offsets) == length(s2.genome.meta_offsets)
        @test length(s1.genome.vocab) == length(s2.genome.vocab)
        @test sort(s1.genome.vocab) == sort(s2.genome.vocab)
    end

    # -------------------------------------------------------------------------
    @testset "add_features! - ID sanitization" begin
        gff = tempname() * ".gff3"
        open(gff, "w") do io
            write(io, "##gff-version 3\n")
            write(
                io,
                "chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\tID=gene:transcript:GENE0001;gene_biotype=protein_coding\n",
            )
        end

        @testset "default sanitization (on)" begin
            s = Species("Test species")
            add_features!(gff, s.genome)

            meta = get_metadata(s.genome, UInt32(1))
            @test length(meta) == 3
            @test meta[1] == "GENE0001"
            @test meta[2] == "RefSeq"
            @test meta[3] == "protein_coding"
        end

        @testset "sanitization disabled" begin
            s = Species("Test species")
            add_features!(gff, s.genome; sanitize_ids = false)

            meta = get_metadata(s.genome, UInt32(1))
            @test length(meta) == 3
            @test meta[1] == "gene:transcript:GENE0001"
            @test meta[2] == "RefSeq"
            @test meta[3] == "protein_coding"
        end
    end

    # -------------------------------------------------------------------------
    @testset "get_metadata – all overloads" begin
        s = Species("C. elegans")
        add_features!(GFF_SINGLE, s.genome)
        scaffold = s.genome.scaffolds["NC_003280.10"]
        n_features = 78_407

        @testset "get_metadata(genome, UInt32) – valid indices" begin
            meta1 = get_metadata(s.genome, UInt32(1))
            @test meta1 isa Vector{String}
            @test length(meta1) == 3   # id, source, biotype

            # Known record verified in the sibling testset above
            @test get_metadata(s.genome, UInt32(2)) ==
                  ["gene-CELE_2L52.2", "RefSeq", "ncRNA"]
        end

        @testset "get_metadata(genome, UInt32) – out-of-range indices" begin
            # length(meta_offsets) == n_features + 1; first invalid index is n_features + 1
            @test get_metadata(s.genome, UInt32(n_features + 1)) == String[]
            @test get_metadata(s.genome, UInt32(n_features + 100)) == String[]
        end

        @testset "get_metadata(genome, IntervalTreeM64)" begin
            results = get_metadata(s.genome, scaffold.features)
            @test results isa Vector{Vector{String}}
            @test length(results) == n_features
            @test all(r -> length(r) == 3, results)
        end

        @testset "get_metadata(genome, Scaffold) – delegates to IntervalTreeM64" begin
            @test get_metadata(s.genome, scaffold) ==
                  get_metadata(s.genome, scaffold.features)
        end

        @testset "get_metadata(genome) – full genome" begin
            all_meta = get_metadata(s.genome)
            @test all_meta isa Dict{String,Vector{Vector{String}}}
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
            @test all(
                iv -> BioinfoTools2.Reference.parse_so_term(iv.value) == gene_code,
                tree,
            )
        end

        @testset "get_feature(scaffold, Symbol) – unknown term" begin
            @test isnothing(get_feature(scaffold, :not_a_real_so_term))
        end

        @testset "get_feature(scaffold, AbstractString) – delegates to Symbol" begin
            @test length(get_feature(scaffold, "gene")) ==
                  length(get_feature(scaffold, :gene))
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

    # -------------------------------------------------------------------------
    @testset "get_so_terms" begin
        @testset "empty genome" begin
            s = Species("C. elegans")
            @test get_so_terms(s.genome) == Symbol[]
        end

        @testset "single-scaffold genome" begin
            s = Species("C. elegans")
            add_features!(GFF_SINGLE, s.genome)

            used_terms = get_so_terms(s.genome)
            @test used_terms isa Vector{Symbol}
            @test !isempty(used_terms)
            @test :gene in used_terms
            @test issorted(used_terms)
        end
    end

    # -------------------------------------------------------------------------
    @testset "getindex – FeatureRecord lookups" begin
        s = Species("C. elegans")
        add_features!(GFF_SINGLE, s.genome)
        scaffold = s.genome.scaffolds["NC_003280.10"]

        # Grab a real gene interval, its metadata index and its ID.
        gene_iv = first(get_feature(scaffold, :gene))
        gene_idx = BioinfoTools2.Reference.parse_index(gene_iv.value)
        gene_id = BioinfoTools2.Reference.get_metadata_id(s.genome, gene_idx)

        @testset "getindex(genome, String) – match / miss" begin
            rec = s.genome[gene_id]
            @test rec isa FeatureRecord
            @test rec.id == gene_id
            @test rec.chromosome == "NC_003280.10"
            @test isnothing(s.genome["NO_SUCH_FEATURE_ID"])
        end

        @testset "getindex(genome, Vector{String})" begin
            recs = s.genome[[gene_id, "NO_SUCH_FEATURE_ID"]]
            @test recs isa Vector{FeatureRecord}
            @test any(r -> r.id == gene_id, recs)
            @test isempty(s.genome[["nope_1", "nope_2"]])
        end

        @testset "getindex(genome, UInt32) – metadata index" begin
            rec = s.genome[gene_idx]
            @test rec isa FeatureRecord
            @test rec.id == gene_id
            @test rec.feature_type == :gene
            @test rec.chromosome == "NC_003280.10"
            @test rec.start_pos == gene_iv.first
            @test rec.end_pos == gene_iv.last

            # An index no feature uses returns nothing.
            @test isnothing(s.genome[typemax(UInt32)])
        end
    end
end
