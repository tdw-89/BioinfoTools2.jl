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

"""Given a full attribute `Dict`, the single value of `key`, or `"NA"`."""
function dict_attribute(attributes::Dict, key::String)
    values = get(attributes, key, String[])
    return length(values) == 1 ? only(values) : "NA"
end

"""
Given a GFF3 file, read its SOFA-typed records the straightforward way — every
attribute decoded into a `Dict`, every type through `convert_so_term`. Returns
`(scaffold, start, end, so_term, strand, metadata)` tuples in file order, which
is what [`loaded_features`](@ref) must reproduce.
"""
function expected_features(path::AbstractString)
    features = Tuple[]
    open(path) do handle
        for record in GFF3.Reader(Reference.gzip_stream(handle, path))
            so_term = convert_so_term(GFF3.featuretype(record))
            isnothing(so_term) && continue
            attributes = Dict(GFF3.attributes(record))
            metadata = [
                sanitize_id(dict_attribute(attributes, "ID")),
                GFF3.hassource(record) ? GFF3.source(record) : "NA",
                dict_attribute(attributes, "gene_biotype"),
            ]
            feature = (
                GFF3.seqname(record),
                UInt32(GFF3.seqstart(record)),
                UInt32(GFF3.seqend(record)),
                so_term,
                GFF3.strand(record),
                metadata,
            )
            push!(features, feature)
        end
    end
    return features
end

"""
Given a loaded genome, list its features by metadata index, shaped like
[`expected_features`](@ref). Returns that vector; an index no interval carries
stays `nothing`.
"""
function loaded_features(genome::Genome)
    features = Vector{Any}(nothing, length(genome.meta_offsets) - 1)
    for (name, scaffold) in genome.scaffolds, interval in scaffold.features
        code = interval.value
        features[parse_index(code)] = (
            name,
            interval.first,
            interval.last,
            parse_so_term(code),
            parse_strand(code),
            get_metadata(genome, parse_index(code)),
        )
    end
    return features
end

"""Given two genomes, compare their metadata stores and every scaffold's intervals."""
function same_genome(left::Genome, right::Genome)
    left.vocab == right.vocab &&
        left.meta_offsets == right.meta_offsets &&
        left.meta_blob == right.meta_blob &&
        keys(left.scaffolds) == keys(right.scaffolds) &&
        all(
            collect(left.scaffolds[name].features) ==
            collect(right.scaffolds[name].features) for name in keys(left.scaffolds)
        )
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
        @test sanitize_id("gene:WBGene00000001") == "WBGene00000001"
        @test sanitize_id("gene:transcript:TX0001") == "TX0001"
        @test sanitize_id("mRNA:RNA:CDS:Protein:ABC123") == "ABC123"
        @test sanitize_id("gene-family:ABC123") == "gene-family:ABC123"
        @test sanitize_id("ABC123") == "ABC123"
    end

    # -------------------------------------------------------------------------
    @testset "parse_record" begin
        @testset "known SO term with sanitized ID and full metadata" begin
            record = read_test_gff_record(
                "chr1\tRefSeq\tgene\t5\t20\t.\t+\t.\tID=gene:transcript:GENE0001;gene_biotype=protein_coding",
            )

            result = parse_record(record, UInt32(7))

            @test result isa ParseResult
            @test result.scaffold_id == "chr1"
            @test result.start_pos == UInt32(5)
            @test result.end_pos == UInt32(20)
            @test parse_index(result.code) == UInt32(7)
            @test parse_strand(result.code) == GFF3.GenomicFeatures.STRAND_POS
            @test parse_so_term(result.code) == convert_so_term("gene")
            @test result.id == "GENE0001"
            @test result.source == "RefSeq"
            @test result.biotype == "protein_coding"
        end

        @testset "sanitization disabled preserves prefixed ID" begin
            record = read_test_gff_record(
                "chr1\tRefSeq\tmRNA\t10\t30\t.\t-\t.\tID=transcript:TX0001;gene_biotype=ncRNA",
            )

            result = parse_record(record, UInt32(8); sanitize_ids = false)

            @test result.id == "transcript:TX0001"
            @test parse_strand(result.code) == GFF3.GenomicFeatures.STRAND_NEG
        end

        @testset "missing optional fields fall back to NA" begin
            record = read_test_gff_record("chr2\t.\tgene\t1\t9\t.\t.\t.\tName=no_id")

            result = parse_record(record, UInt32(9))

            @test result.id == "NA"
            @test result.source == "NA"
            @test result.biotype == "NA"
            @test parse_strand(result.code) == GFF3.GenomicFeatures.STRAND_BOTH
        end

        @testset "multi-valued ID and biotype fall back to NA" begin
            record = read_test_gff_record(
                "chr3\tRefSeq\tgene\t2\t8\t.\t+\t.\tID=gene:ONE,gene:TWO;gene_biotype=type1,type2",
            )

            result = parse_record(record, UInt32(10))

            @test result.id == "NA"
            @test result.source == "RefSeq"
            @test result.biotype == "NA"
        end

        @testset "unknown SO term is skipped" begin
            record = read_test_gff_record(
                "chr4\tRefSeq\tnot_a_real_so_term\t3\t6\t.\t+\t.\tID=gene:SKIPME;gene_biotype=protein_coding",
            )

            @test isnothing(parse_record(record, UInt32(11)))
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
        @test all(id -> haskey(s.genome.scaffolds, id), expected_scaffolds)

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
    @testset "attribute_or_na – one attribute, read as a full decode would" begin
        keys_read = ["ID", "gene_biotype", "Name", "Parent", "no_such_key"]
        matches = Bool[]
        open(GFF_SINGLE) do handle
            for record in GFF3.Reader(Reference.gzip_stream(handle, GFF_SINGLE))
                attributes = Dict(GFF3.attributes(record))
                push!(
                    matches,
                    all(
                        Reference.attribute_or_na(record, key) ==
                        dict_attribute(attributes, key) for key in keys_read
                    ),
                )
            end
        end
        @test length(matches) == 78_559   # every record, SOFA-typed or not
        @test all(matches)

        # Percent-decoding after the comma split, a repeated key read last (as a
        # `Dict` keeps it), and a multi-valued attribute.
        record = read_test_gff_record(
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=one;Note=a%2Cb,c;gene_biotype=x%3By;ID=two",
        )
        @test Reference.attribute_or_na(record, "ID") == "two"
        @test Reference.attribute_or_na(record, "gene_biotype") == "x;y"
        @test Reference.attribute_or_na(record, "Note") == "NA"
    end

    # -------------------------------------------------------------------------
    @testset "add_features! – feature_types and biotypes" begin
        everything = expected_features(GFF_SINGLE)
        gene, mrna = convert_so_term("gene"), convert_so_term("mRNA")
        function filtered_load(; options...)
            genome = Species("C. elegans").genome
            add_features!(GFF_SINGLE, genome; options...)
            return loaded_features(genome)
        end

        # Unfiltered, the loader reproduces the straightforward reading record for
        # record, metadata indices staying dense in file order.
        @test filtered_load() == everything
        @test filtered_load(feature_types = :gene) == filter(f -> f[4] == gene, everything)
        @test filtered_load(feature_types = [:gene, "mRNA"]) ==
              filter(f -> f[4] in (gene, mrna), everything)

        protein_coding =
            filter(f -> f[4] == gene && f[6][3] == "protein_coding", everything)
        @test length(protein_coding) == 3_446
        @test filtered_load(feature_types = :gene, biotypes = "protein_coding") ==
              protein_coding
        # A record without gene_biotype reads as "NA", so only listing "NA" keeps it.
        @test filtered_load(biotypes = [:protein_coding, "NA"]) ==
              filter(f -> f[6][3] in ("protein_coding", "NA"), everything)
        @test isempty(filtered_load(feature_types = Symbol[]))

        species = Species("C. elegans")
        add_features!(GFF_SINGLE, species; feature_types = :gene, sanitize_ids = false)
        @test get_so_terms(species.genome) == [:gene]
        @test length(loaded_features(species.genome)) == 4_893

        # SOFA has no lncRNA term (NCBI writes `lncRNA` and `lnc_RNA`).
        @test_throws ArgumentError add_features!(
            GFF_SINGLE,
            Species("C. elegans").genome;
            feature_types = [:gene, :lncRNA],
        )
    end

    # -------------------------------------------------------------------------
    @testset "add_features! – successive files share one metadata store" begin
        # Alternating scaffolds and sources, one of them missing.
        extra = """
            chrNew\tsrcA\tgene\t5\t50\t.\t-\t.\tID=gene:NEW1;gene_biotype=protein_coding
            NC_003280.10\t.\tgene\t100\t200\t.\t+\t.\tID=NEW2;gene_biotype=lncRNA
            chrNew\tsrcB\tmRNA\t5\t40\t.\t-\t.\tID=rna-NEW1;Parent=gene:NEW1
            """
        second = tempname() * ".gff3"
        combined = tempname() * ".gff3"
        try
            write(second, "##gff-version 3\n" * extra)
            first_text = open(GFF_SINGLE) do handle
                read(Reference.gzip_stream(handle, GFF_SINGLE), String)
            end
            write(combined, first_text * extra)

            successive = Species("C. elegans").genome
            add_features!(GFF_SINGLE, successive)
            add_features!(second, successive)
            at_once = Species("C. elegans").genome
            add_features!(combined, at_once)

            @test same_genome(successive, at_once)
            @test successive["NEW2"].metadata == ["NEW2", "NA", "lncRNA"]
            @test successive["rna-NEW1"].chromosome == "chrNew"
        finally
            rm(second; force = true)
            rm(combined; force = true)
        end
    end

    # -------------------------------------------------------------------------
    @testset "add_features! – directives and comments after the last record" begin
        records = """
            ##gff-version 3
            chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=geneA;gene_biotype=protein_coding
            ##sequence-region chr2 1 100
            chr2\tsrc\tgene\t5\t20\t.\t-\t.\tID=geneB;gene_biotype=protein_coding
            """
        # NCBI files close with `###`; GFF3.jl used to throw EOFError on it.
        endings = ["###\n", "# a closing comment\n", "##FASTA\n>chr1\nACGTACGTA\n"]
        loaded = map(endings) do ending
            gff = tempname() * ".gff3"
            try
                write(gff, records * ending)
                genome = Species("Test species").genome
                add_features!(gff, genome)
                [record.id for record in genome[["geneA", "geneB"]]]
            finally
                rm(gff; force = true)
            end
        end
        @test all(ids -> sort(ids) == ["geneA", "geneB"], loaded)

        # A record cut off mid-line is still an error, not a silent partial load.
        gff = tempname() * ".gff3"
        try
            write(gff, records * "chr2\tsrc\tge")
            @test_throws EOFError add_features!(gff, Species("Test species").genome)
        finally
            rm(gff; force = true)
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
            gene_code = parse_so_term(first(tree).value)
            @test all(iv -> parse_so_term(iv.value) == gene_code, tree)
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
            @test all(
                length(by_symbol[name]) == length(by_string[name]) for
                name in keys(by_symbol)
            )
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
        gene_idx = parse_index(gene_iv.value)
        gene_id = get_metadata_id(s.genome, gene_idx)

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

        @testset "getindex(genome, Vector{UInt32}) – batched metadata indices" begin
            # Every gene's index resolved in one walk must match what the scalar
            # method returns index by index.
            gene_indices =
                UInt32[parse_index(iv.value) for iv in get_feature(scaffold, :gene)]
            @test length(gene_indices) > 1

            batched = s.genome[gene_indices]
            @test batched isa Dict{UInt32,FeatureRecord}
            @test length(batched) == length(gene_indices)

            matches = map(gene_indices) do idx
                scalar = s.genome[idx]
                record = batched[idx]
                record.id == scalar.id &&
                    record.feature_type == scalar.feature_type &&
                    record.chromosome == scalar.chromosome &&
                    record.start_pos == scalar.start_pos &&
                    record.end_pos == scalar.end_pos &&
                    record.code == scalar.code
            end
            @test all(matches)

            # Unused indices are simply absent; an empty request walks nothing.
            mixed = s.genome[UInt32[gene_idx, typemax(UInt32)]]
            @test collect(keys(mixed)) == [gene_idx]
            @test isempty(s.genome[UInt32[]])
            @test isempty(s.genome[UInt32[typemax(UInt32)]])
        end
    end

    @testset "feature_id" begin
        species = Species("C. elegans")
        add_features!(GFF_SINGLE, species.genome)
        scaffold = first(values(species.genome.scaffolds))
        intervals = collect(scaffold.features)

        ids = [feature_id(species.genome, interval) for interval in intervals[1:20]]
        @test all(!isnothing, ids)
        # It resolves exactly what the record carries, without rebuilding one.
        records = [
            feature_record(species.genome, scaffold.name, interval) for
            interval in intervals[1:20]
        ]
        @test ids == [record.id for record in records]

        # An interval whose metadata index points past the store has no ID.
        stray = Reference.IntervalValue(UInt32(1), UInt32(2), UInt64(typemax(UInt32)))
        @test feature_id(species.genome, stray) === nothing
    end
end
