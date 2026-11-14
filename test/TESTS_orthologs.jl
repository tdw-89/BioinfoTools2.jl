using BioinfoTools2.Reference
using BioinfoTools2.Orthologs
using DataFrames
using Test

const OR = BioinfoTools2.Homologs.Orthologs

# Four species: `out` has no genome, and `big` is a prefix of `big_ape`, so leaves
# must be split on the longest species name. Leaves in preorder: out, ape, big, big_ape.
const OR_SPECIES_TREE = "(out:0.5,(ape:0.2,(big:0.1,big_ape:0.1)N2:0.05)N1:0.3)N0;"

# OG0000001 holds a duplication (n1), OG0000002 a polytomy root whose children
# overlap only in part, and OG0000005 a protein (A9.1) that no GFF codes for.
const OR_GENE_TREES = [
    "OG0000000: (out_O1.1:0.4,(ape_A1.1:0.2,(big_B1.1:0.1,big_ape_C1.1:0.1)n2:0.05)n1:0.3)n0;",
    "OG0000001: (out_O2.1:0.5,((ape_A2.1:0.1,big_B2.1:0.1)n2:0.1,(ape_A3.1:0.1,big_B3.1:0.2)n3:0.1)n1:0.2)n0;",
    "OG0000002: (ape_A4.1:0.1,big_B4.1:0.1,(big_B5.1:0.1,big_ape_C2.1:0.1)n1:0.1,out_O3.1:0.3)n0;",
    "OG0000005: (ape_A9.1:0.1,big_ape_C3.1:0.1)n0;",
]

# Columns in another order than the tree; OG0000003 (three genes) and OG0000004
# (two) have no tree, and OG0000004's row stops short of its empty cells.
const OR_ORTHOGROUPS = [
    "Orthogroup\tape\tbig\tbig_ape\tout",
    "OG0000000\tA1.1\tB1.1\tC1.1\tO1.1",
    "OG0000001\tA2.1, A3.1\tB2.1, B3.1\t\tO2.1",
    "OG0000002\tA4.1\tB4.1, B5.1\tC2.1\tO3.1",
    "OG0000003\tA5.1\tB6.1\t\tO4.1",
    "OG0000004\tA6.1\tB7.1",
    "OG0000005\tA9.1\t\tC3.1\t",
]

"""Given NCBI-style IDs, one gene's GFF3 lines: gene → mRNA → exon/CDS, or a CDS parented by the gene when `direct`."""
function ncbi_gene(gene_id, protein, start; direct = false, symbol = nothing)
    stop, rna = start + 99, "rna-$protein"
    gene_attributes = isnothing(symbol) ? "ID=$gene_id" : "ID=$gene_id;gene=$symbol"
    lines = [
        "chr1\tRefSeq\tgene\t$start\t$stop\t.\t+\t.\t$gene_attributes;gene_biotype=protein_coding",
    ]
    if !direct
        push!(lines, "chr1\tRefSeq\tmRNA\t$start\t$stop\t.\t+\t.\tID=$rna;Parent=$gene_id")
        push!(
            lines,
            "chr1\tRefSeq\texon\t$start\t$stop\t.\t+\t.\tID=exon-$protein-1;Parent=$rna",
        )
    end
    parent = direct ? gene_id : rna
    push!(
        lines,
        "chr1\tRefSeq\tCDS\t$(start+10)\t$(stop-10)\t.\t+\t0\tID=cds-$protein;Parent=$parent;protein_id=$protein",
    )
    return lines
end

"""Given Ensembl-style (prefixed) IDs, one gene's GFF3 lines: gene → transcript → CDS."""
ensembl_gene(name, start) = [
    "chr2\tensembl\tgene\t$start\t$(start+99)\t.\t-\t.\tID=gene:$name;biotype=protein_coding",
    "chr2\tensembl\tmRNA\t$start\t$(start+99)\t.\t-\t.\tID=transcript:T$name;Parent=gene:$name",
    "chr2\tensembl\tCDS\t$(start+10)\t$(start+89)\t.\t-\t0\tID=CDS:$name.1;Parent=transcript:T$name;protein_id=$name.1",
]

"""Given GFF3 record lines, write them under a header and NCBI's closing `###`. Returns the path."""
function write_gff(directory, name, lines)
    path = joinpath(directory, "$name.gff3")
    write(path, "##gff-version 3\n" * join(lines, '\n') * "\n###\n")
    return path
end

"""
Given the gene-tree and orthogroup lines, write an OrthoFinder-shaped results
directory and one GFF3 per genome species. Returns `(directory, gffs)`.
"""
function clade_fixture(; gene_trees = OR_GENE_TREES, orthogroups = OR_ORTHOGROUPS)
    directory = mktempdir()
    for folder in ("Species_Tree", "Resolved_Gene_Trees", "Orthogroups")
        mkpath(joinpath(directory, folder))
    end
    write(
        joinpath(directory, "Species_Tree", "SpeciesTree_rooted_node_labels.txt"),
        OR_SPECIES_TREE,
    )
    write(
        joinpath(directory, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt"),
        join(gene_trees, '\n') * '\n',
    )
    write(
        joinpath(directory, "Orthogroups", "Orthogroups.tsv"),
        join(orthogroups, '\n') * '\n',
    )

    # ape: A3's gene ID is suffixed (`gene-SHARED-2`), A5's CDS hangs off the gene.
    ape_lines = vcat(
        [ncbi_gene("gene-A$i", "A$i.1", 1000i) for i in (1, 2, 4, 6)]...,
        ncbi_gene("gene-SHARED-2", "A3.1", 3000; symbol = "SHARED"),
        ncbi_gene("gene-A5", "A5.1", 5000; direct = true),
    )
    gffs = Dict(
        "ape" => write_gff(directory, "ape", ape_lines),
        "big" => write_gff(
            directory,
            "big",
            vcat([ensembl_gene("B$i", 1000i) for i = 1:7]...),
        ),
        "big_ape" => write_gff(
            directory,
            "big_ape",
            vcat([ncbi_gene("gene-C$i", "C$i.1", 1000i) for i = 1:3]...),
        ),
    )
    return directory, gffs
end

"""Given a clade's ortholog table, its rows as `(orthogroup, protein_a, protein_b)` triples."""
pair_set(table) = Set(zip(table.orthogroup, table.protein_a, table.protein_b))

@testset "Orthologs" begin
    directory, gffs = clade_fixture()
    clade =
        @test_logs (:warn, r"^ape: 1 of 7 gene-tree proteins have no gene \(1 not in") match_mode =
            :any load_clade(directory, gffs)
    forest = clade.gene_trees

    # -------------------------------------------------------------------------
    @testset "_newick_tokens lengths" begin
        tokens = BioinfoTools2.Homologs.Paralogs._newick_tokens
        @test tokens("(A:1,B:2.5e-3)n0:4;"; lengths = true) ==
              ['(', "A", 1.0, "B", 0.0025, ')', "n0", 4.0]
        # The default stays what Paralogs' depth and leaf readers expect.
        @test tokens("(A:1,B:2.5e-3)n0:4;") == ['(', "A", "B", ')', "n0"]
        @test tokens("(A:1,B:2)"; lengths = true) == ['(', "A", 1.0, "B", 2.0, ')']
    end

    # -------------------------------------------------------------------------
    @testset "SpeciesTree" begin
        tree = clade.species_tree
        @test species_names(tree) == ["out", "ape", "big", "big_ape"]
        @test tree.labels == ["N0", "out", "N1", "ape", "N2", "big", "big_ape"]
        # Depths number as lca_depth does (Paralogs' reading of the same tree).
        @test Dict(zip(tree.labels, tree.depth)) ==
              BioinfoTools2.Homologs.Paralogs._newick_depths(OR_SPECIES_TREE)
        @test tree.species == UInt64[0b1111, 0b0001, 0b1110, 0b0010, 0b1100, 0b0100, 0b1000]
        @test tree.branch_length[2] == 0.5f0 && isnan(tree.branch_length[1])

        # SpeciesTree_rooted.txt carries support values; its nodes get the same names.
        supported = SpeciesTree("(out:0.5,(ape:0.2,(big:0.1,big_ape:0.1)1:0.05)0.98:0.3);")
        @test supported.labels == tree.labels
        @test_throws ArgumentError SpeciesTree("(out,(ape,ape)N1)N0;")
        @test_throws ArgumentError SpeciesTree("((out,ape)N1;")
        @test sprint(show, tree) == "SpeciesTree(out, ape, big, big_ape)"
    end

    # -------------------------------------------------------------------------
    @testset "protein_gene_map" begin
        @test protein_gene_map(gffs["ape"]) == Dict(
            "A1.1" => "gene-A1",
            "A2.1" => "gene-A2",
            "A3.1" => "gene-SHARED-2",   # not "gene-" * symbol
            "A4.1" => "gene-A4",
            "A5.1" => "gene-A5",         # CDS parented by the gene itself
            "A6.1" => "gene-A6",
        )
        # IDs are sanitized as add_features! sanitizes them, unless told not to.
        @test protein_gene_map(gffs["big"]) == Dict("B$i.1" => "B$i" for i = 1:7)
        @test protein_gene_map(gffs["big"]; sanitize_ids = false) ==
              Dict("B$i.1" => "gene:B$i" for i = 1:7)
    end

    # -------------------------------------------------------------------------
    @testset "Clade layout" begin
        @test species_names(clade) == ["out", "ape", "big", "big_ape"]
        @test isempty(clade.species[1].genome.scaffolds)   # `out` has no GFF
        @test orthogroups(clade) == ["OG000000$i" for i = 0:5]
        @test forest.source == [
            OR.TREE_RESOLVED,
            OR.TREE_RESOLVED,
            OR.TREE_RESOLVED,
            OR.TREE_STAR,
            OR.TREE_PAIR,
            OR.TREE_RESOLVED,
        ]

        # Preorder invariants: parents come first, a subtree is the range up to
        # its subtree_end, each child's range nests inside its parent's, and an
        # internal mask is exactly the union of its children's.
        invariants = Bool[]
        for tree in eachindex(forest.orthogroups), node in OR._tree_nodes(forest, tree)
            tree_range = OR._tree_nodes(forest, tree)
            parent = forest.parent[node]
            kids = collect(OR._Children(forest.subtree_end, node))
            push!(
                invariants,
                (parent == 0) == (node == first(tree_range)) &&
                    (parent == 0 || first(tree_range) <= parent < node) &&
                    node <= forest.subtree_end[node] <= last(tree_range) &&
                    all(forest.parent[kid] == node for kid in kids) &&
                    (
                        isempty(kids) ? count_ones(forest.species[node]) == 1 :
                        forest.species[node] == reduce(|, forest.species[kids]) &&
                        forest.subtree_end[node] == forest.subtree_end[last(kids)]
                    ),
            )
        end
        @test length(invariants) == length(forest.parent) == 33
        @test all(invariants)

        # Resolved trees write back exactly; a star gets OrthoFinder's `n0` root.
        trees = Dict(split(line, ": "; limit = 2) for line in OR_GENE_TREES)
        @test all(newick(gene_tree(clade, og)) == text for (og, text) in trees)
        @test newick(gene_tree(clade, "OG0000003")) == "(ape_A5.1,big_B6.1,out_O4.1)n0;"
        @test sprint(show, clade) == "Clade(out, ape, big, big_ape; 6 orthogroups)"
        @test sprint(show, forest) == "GeneForest(6 trees, 33 nodes)"
    end

    # -------------------------------------------------------------------------
    @testset "orthogroup" begin
        @test orthogroup(clade, "ape", "gene-SHARED-2") == "OG0000001"
        @test orthogroup(clade, "ape", "A1.1") == "OG0000000"       # by protein
        @test orthogroup(clade, "out", "O3.1") == "OG0000002"       # genome-less
        @test orthogroup(clade, "ape", "A9.1") == "OG0000005"       # unplaced protein
        @test orthogroup(clade, "ape", "gene-A6") == "OG0000004"    # two-gene orthogroup
        @test isnothing(orthogroup(clade, "ape", "gene-NOPE"))
        @test_throws ArgumentError orthogroup(clade, "cat", "gene-A1")
    end

    # -------------------------------------------------------------------------
    @testset "orthologs" begin
        found = orthologs(clade, "ape", "gene-A2")
        @test found.species == ["big", "out"]
        @test isequal(found.gene, ["B2", missing])
        @test found.protein == ["B2.1", "O2.1"]
        # O2 is orthologous to both A2 and A3, so that pair is not one-to-one.
        @test found.one_to_one == [true, false]

        # A genome-less species is queried by protein; the duplication at n1 puts
        # all four genes below it.
        @test orthologs(clade, "out", "O2.1").protein == ["A2.1", "B2.1", "A3.1", "B3.1"]

        # The polytomy decides per child pair: B4 shares `big` with (B5, C2), and
        # A4 and O3 each have both B4 and B5 as orthologs in `big`.
        polytomy = orthologs(clade, "big", "B4")
        @test polytomy.protein == ["A4.1", "O3.1"]
        @test polytomy.one_to_one == [false, false]
        @test orthologs(clade, "big_ape", "gene-C2").protein == ["B5.1", "A4.1", "O3.1"]
        @test orthologs(clade, "big_ape", "gene-C2"; species_wanted = "ape").protein ==
              ["A4.1"]
        @test orthologs(clade, "big_ape", "gene-C2"; species_wanted = ["ape", "out"]).protein ==
              ["A4.1", "O3.1"]

        # A star calls every cross-species pair; a two-gene orthogroup calls none.
        star = orthologs(clade, "ape", "gene-A5")
        @test star.protein == ["B6.1", "O4.1"] && all(star.one_to_one)
        @test nrow(orthologs(clade, "ape", "gene-A6")) == 0
        @test nrow(orthologs(clade, "ape", "gene-NOPE")) == 0
        @test names(orthologs(clade, "ape", "gene-NOPE")) ==
              ["species", "gene", "protein", "one_to_one"]
        @test_throws ArgumentError orthologs(
            clade,
            "ape",
            "gene-A1";
            species_wanted = "cat",
        )
    end

    # -------------------------------------------------------------------------
    @testset "ortholog_pairs" begin
        pairs = ortholog_pairs(clade, "ape", "big")
        @test pair_set(pairs) == Set([
            ("OG0000000", "A1.1", "B1.1"),
            ("OG0000001", "A2.1", "B2.1"),
            ("OG0000001", "A3.1", "B3.1"),
            ("OG0000002", "A4.1", "B4.1"),
            ("OG0000002", "A4.1", "B5.1"),
            ("OG0000003", "A5.1", "B6.1"),
        ])
        @test nrow(pairs) == 6
        @test Set(pairs.protein_a[.!pairs.one_to_one]) == Set(["A4.1"])
        @test pairs.gene_a[pairs.protein_a .== "A3.1"] == ["gene-SHARED-2"]

        # Every pair the table lists is one orthologs() finds, gene by gene, for
        # every species pair.
        consistent = Bool[]
        for a in species_names(clade), b in species_names(clade)
            a == b && continue
            table = ortholog_pairs(clade, a, b)
            by_gene = Set(
                (orthogroup(clade, a, protein), protein, other) for
                protein in clade.proteins[OR._species_position(clade, a)].keys for
                other in orthologs(clade, a, protein; species_wanted = b).protein
            )
            push!(consistent, pair_set(table) == by_gene)
        end
        @test length(consistent) == 12 && all(consistent)

        # Mirrored, and with the genome-less species' genes missing.
        mirrored = ortholog_pairs(clade, "big", "ape")
        @test pair_set(mirrored) == Set((og, b, a) for (og, a, b) in pair_set(pairs))
        @test all(ismissing, ortholog_pairs(clade, "ape", "out").gene_b)
        @test_throws ArgumentError ortholog_pairs(clade, "ape", "ape")
    end

    # -------------------------------------------------------------------------
    @testset "common_ancestor" begin
        paralogs = common_ancestor(clade, "ape" => "gene-A2", "ape" => "gene-SHARED-2")
        @test paralogs.event == :duplication
        @test (paralogs.species_node, paralogs.depth) == ("N1", 2)
        @test node_label(gene_tree(clade, "OG0000001"), paralogs.node) == "n1"

        @test common_ancestor(clade, "ape" => "gene-A2", "big" => "B2").event == :speciation
        @test common_ancestor(clade, "ape" => "gene-A2", "big" => "B3").event ==
              :duplication

        # At the polytomy, two pairs meeting at one node differ in their event.
        overlapping = common_ancestor(clade, "big" => "B4", "big_ape" => "gene-C2")
        disjoint = common_ancestor(clade, "ape" => "gene-A4", "big_ape" => "gene-C2")
        @test overlapping.node == disjoint.node
        @test (overlapping.event, disjoint.event) == (:duplication, :speciation)
        @test (overlapping.species_node, overlapping.depth) == ("N0", 1)

        @test common_ancestor(clade, "out" => "O1.1", "big_ape" => "gene-C1").event ==
              :speciation
        @test common_ancestor(clade, "ape" => "gene-A6", "big" => "B7").event == :unresolved
        @test isnothing(common_ancestor(clade, "ape" => "gene-A1", "big" => "B2"))
        @test isnothing(common_ancestor(clade, "ape" => "gene-NOPE", "big" => "B2"))
        @test_throws ArgumentError common_ancestor(
            clade,
            "ape" => "gene-A1",
            "ape" => "A1.1",
        )
    end

    # -------------------------------------------------------------------------
    @testset "gene_tree access" begin
        tree = gene_tree(clade, "OG0000001")
        top = root(tree)
        @test length(tree_nodes(tree)) == 9 && first(tree_nodes(tree)) == top
        @test isnothing(parent(tree, top)) && isnan(branch_length(tree, top))
        kids = children(tree, top)
        @test node_label.(Ref(tree), kids) == ["out_O2.1", "n1"]
        @test parent(tree, kids[2]) == top && branch_length(tree, kids[2]) == 0.2f0
        @test is_leaf(tree, kids[1]) && !is_leaf(tree, kids[2])
        @test node_label.(Ref(tree), leaves(tree)) ==
              ["out_O2.1", "ape_A2.1", "big_B2.1", "ape_A3.1", "big_B3.1"]
        @test length(leaves(tree, kids[2])) == 4

        shared = leaves(tree)[4]
        @test (
            leaf_species(tree, shared),
            leaf_gene(tree, shared),
            leaf_protein(tree, shared),
        ) == ("ape", "gene-SHARED-2", "A3.1")
        @test ismissing(leaf_gene(tree, kids[1]))
        unplaced = gene_tree(clade, "OG0000005")
        @test ismissing(leaf_gene(unplaced, first(leaves(unplaced))))
        @test sprint(show, tree) == "GeneTree(\"OG0000001\", 5 leaves)"

        @test_throws ArgumentError gene_tree(clade, "OG9999999")
        @test_throws ArgumentError children(tree, root(unplaced))   # another tree's node
    end

    # -------------------------------------------------------------------------
    @testset "Clade from loaded species" begin
        loaded = [Species(name) for name in ("ape", "big")]
        add_features!(gffs["ape"], loaded[1]; feature_types = :gene)
        add_features!(gffs["big"], loaded[2]; feature_types = :gene)
        maps = Dict(name => protein_gene_map(gffs[name]) for name in ("ape", "big"))
        rebuilt = @test_logs (:warn, r"^ape: 1 of 7") Clade(
            loaded,
            directory;
            protein_maps = maps,
        )
        @test pair_set(ortholog_pairs(rebuilt, "ape", "big")) ==
              pair_set(ortholog_pairs(clade, "ape", "big"))
        # big_ape was not handed in, so it is genome-less here.
        @test all(ismissing, ortholog_pairs(rebuilt, "ape", "big_ape").gene_b)

        @test_logs (:warn, r"big has a genome but no protein map") match_mode = :any Clade(
            loaded,
            directory;
            protein_maps = Dict("ape" => maps["ape"]),
        )
        # Without Orthogroups.tsv, the tree-less orthogroups are simply absent.
        resolved_only = @test_logs (:warn, r"^ape: 1 of 5") Clade(
            loaded,
            directory;
            protein_maps = maps,
            orthogroups = nothing,
        )
        @test length(orthogroups(resolved_only)) == 4
    end

    # -------------------------------------------------------------------------
    @testset "malformed input" begin
        # A species the tree lacks throws before any GFF loads.
        @test_throws ArgumentError load_clade(directory, Dict("cat" => gffs["ape"]))
        @test_throws ArgumentError Clade([Species("ape"), Species("ape")], directory)

        stray_leaf =
            clade_fixture(gene_trees = ["OG0000000: (cat_X1.1:0.1,ape_A1.1:0.1)n0;"])[1]
        @test_throws ArgumentError Clade(Species[], stray_leaf)
        misnumbered =
            clade_fixture(gene_trees = ["OG0000000: ((ape_A1.1,big_B1.1)n5,out_O1.1)n0;"])[1]
        @test_throws ArgumentError Clade(Species[], misnumbered)
        unbalanced = clade_fixture(gene_trees = ["OG0000000: ((ape_A1.1,big_B1.1)n1;"])[1]
        @test_throws ArgumentError Clade(Species[], unbalanced)
        unlabelled = clade_fixture(gene_trees = ["OG0000000: (:0.1,ape_A1.1:0.2)n0;"])[1]
        @test_throws ArgumentError Clade(Species[], unlabelled)
        malformed = clade_fixture(gene_trees = ["OG0000000 (ape_A1.1,big_B1.1)n0;"])[1]
        @test_throws ArgumentError Clade(Species[], malformed)
        stray_column =
            clade_fixture(orthogroups = ["Orthogroup\tape\tcat", "OG0000009\tA1.1\tX1.1"])[1]
        @test_throws ArgumentError Clade(Species[], stray_column)
    end

    # -------------------------------------------------------------------------
    @testset "node_event and species_node on gene trees" begin
        duplicated = gene_tree(clade, "OG0000001")
        top, inner = root(duplicated), children(duplicated, root(duplicated))[2]
        @test node_event.(Ref(duplicated), (top, inner, leaves(duplicated)[1])) ==
              (:speciation, :duplication, :leaf)
        @test species_node(duplicated, inner) == "N1" &&
              species_node(duplicated, top) == "N0"
        polytomy = gene_tree(clade, "OG0000002")
        @test node_event(polytomy, root(polytomy)) == :mixed
        pair = gene_tree(clade, "OG0000004")
        @test node_event(pair, root(pair)) == :unresolved
        @test (
            orthogroup(clade, "C2.1"),
            orthogroup(clade, "O3.1"),
            orthogroup(clade, "gene-A1"),
            orthogroup(clade, "nope"),
        ) == ("OG0000002", "OG0000002", "OG0000000", nothing)
    end

    # -------------------------------------------------------------------------
    # OG0000001 as codeml writes it: unrooted, at a trifurcation over n2, n3 and O2.1.
    divergence_dS = "((A2.1:0.3,B2.1:0.1):0.2,(A3.1:0.2,B3.1:0.4):0.6,O2.1:1.4);"
    divergence_dN = "((A2.1:0.15,B2.1:0.05):0.1,(A3.1:0.1,B3.1:0.2):0.3,O2.1:0.7);"

    @testset "divergence_tree rooting and events" begin
        tree = divergence_tree(clade, "OG0000001"; dS = divergence_dS, dN = divergence_dN)
        tips, top = leaves(tree), children(tree, root(tree))
        @test leaf_protein.(Ref(tree), tips) == ["O2.1", "A2.1", "B2.1", "A3.1", "B3.1"]
        # OrthoFinder's lengths below its root are 0.5 (O2.1) and 0.2 (n1): 1.4 splits 1.0 + 0.4.
        @test branch_length.(Ref(tree), top) ≈ [1.0, 0.4] &&
              branch_length(tree, top[1]; rate = :dN) ≈ 0.5
        @test isnothing(parent(tree, root(tree))) && isnan(branch_length(tree, root(tree)))
        @test node_event.(Ref(tree), tree_nodes(tree)) == [
            :speciation,
            :leaf,
            :duplication,
            :speciation,
            :leaf,
            :leaf,
            :speciation,
            :leaf,
            :leaf,
        ]
        @test all(tree.exact) && isempty(tree.dropped)
        # Lengths are stored as Float32, hence the tolerance.
        @test isapprox(patristic(tree, tips[2], tips[5]), 1.5; atol = 1e-6) &&
              isapprox(patristic(tree, tips[2], tips[5]; rate = :dN), 0.75; atol = 1e-6)
        @test isapprox(root_distances(tree)[tips[2]], 0.9; atol = 1e-6)
        @test (species_node(tree, top[2]), species_node(tree, tips[2])) == ("N1", "ape")
        @test (leaf_species(tree, tips[2]), leaf_gene(tree, tips[2])) == ("ape", "gene-A2")
        @test ismissing(leaf_gene(tree, tips[1]))
        @test sprint(show, tree) == "DivergenceTree(\"OG0000001\", 5 leaves)"

        # The same unrooted tree, written rooted on the O2.1 branch, lands identically.
        bifurcating = "(O2.1:0.7,((A2.1:0.3,B2.1:0.1):0.2,(A3.1:0.2,B3.1:0.4):0.6):0.7);"
        rerooted = divergence_tree(clade, "OG0000001"; dS = bifurcating)
        @test rerooted.parent == tree.parent &&
              rerooted.dS[2:end] ≈ tree.dS[2:end] &&
              all(isnan, rerooted.dN)

        # Two genes, no OrthoFinder lengths: the root falls midway, the event unresolved.
        pair = divergence_tree(clade, "OG0000004"; dS = "(A6.1:0.2,B7.1:0.6);")
        @test branch_length.(Ref(pair), children(pair, 1)) ≈ [0.4, 0.4] &&
              node_event(pair, 1) == :unresolved
    end

    # -------------------------------------------------------------------------
    @testset "divergence_tree at a polytomy and with dropped genes" begin
        # codeml resolves OG0000002's four-way root as ((A4, B4), n1, O3).
        resolved = divergence_tree(
            clade,
            "OG0000002";
            dS = "((A4.1:0.1,B4.1:0.2):0.3,(B5.1:0.1,C2.1:0.1):0.2,O3.1:0.5);",
        )
        @test leaf_protein.(Ref(resolved), leaves(resolved)) ==
              ["A4.1", "B4.1", "B5.1", "C2.1", "O3.1"]
        @test node_event.(Ref(resolved), tree_nodes(resolved)) ==
              [:mixed, :speciation, :leaf, :leaf, :speciation, :leaf, :leaf, :leaf]
        @test resolved.exact == [true, false, true, true, true, true, true, true]
        @test species_node(resolved, 2) == "N1"

        # Without B3.1, A3.1 still stands for n3's lineage, so n1 stays a duplication.
        partial = divergence_tree(
            clade,
            "OG0000001";
            dS = "((A2.1:0.3,B2.1:0.1):0.2,A3.1:0.8,O2.1:1.4);",
        )
        lone = last(leaves(partial))
        @test forest.protein[partial.dropped] == ["B3.1"]
        @test node_event(partial, parent(partial, lone)) == :duplication
        @test species_node(partial, lone) == "ape"
        @test prune(partial, ["ape", "big"]).dropped == partial.dropped &&
              isempty(prune(partial, "ape").dropped)
    end

    # -------------------------------------------------------------------------
    @testset "prune" begin
        tree = divergence_tree(clade, "OG0000001"; dS = divergence_dS, dN = divergence_dN)
        primates = prune(tree, ["ape", "big"])
        @test leaf_protein.(Ref(primates), leaves(primates)) ==
              ["A2.1", "B2.1", "A3.1", "B3.1"]
        @test node_event.(Ref(primates), tree_nodes(primates)) ==
              [:duplication, :speciation, :leaf, :leaf, :speciation, :leaf, :leaf]
        # The root lost O2.1, so n1 is the new root and its branch is the stem.
        @test branch_length(primates, 1) ≈ 0.4 &&
              branch_length(primates, 1; rate = :dN) ≈ 0.2
        apes = prune(tree, "ape")
        @test branch_length.(Ref(apes), leaves(apes)) ≈ [0.5, 0.8] &&
              node_event(apes, 1) == :duplication
        @test isapprox(patristic(apes, leaves(apes)...), 1.3; atol = 1e-6)
        @test prune(apes, "ape").parent == apes.parent
        @test_throws ArgumentError prune(tree, "out")
    end

    # -------------------------------------------------------------------------
    @testset "CodemlTrees and malformed codeml trees" begin
        codeml = mktempdir()
        write(joinpath(codeml, "OG0000001_dS.nwk"), divergence_dS)
        write(joinpath(codeml, "OG0000001_dN.nwk"), divergence_dN)
        write(
            joinpath(codeml, "OG0000002_dS.nwk"),
            "((A4.1:0.1,B4.1:0.2):0.3,(B5.1:0.1,C2.1:0.1):0.2,O3.1:0.5);",
        )
        write(joinpath(codeml, "OG0000001_ks.tsv"), "")
        trees = CodemlTrees(codeml)
        @test orthogroups(trees) == ["OG0000001", "OG0000002"] && "OG0000002" in trees
        @test branch_length(divergence_tree(clade, trees, "OG0000001"), 2; rate = :dN) ≈ 0.5
        @test all(isnan, divergence_tree(clade, trees, "OG0000002").dN)
        @test_throws ArgumentError divergence_tree(clade, trees, "OG0000003")
        @test_throws ArgumentError CodemlTrees(mktempdir())

        conflicting = "((A2.1:0.3,A3.1:0.1):0.2,(B2.1:0.2,B3.1:0.4):0.6,O2.1:1.4);"
        @test_throws ArgumentError divergence_tree(clade, "OG0000001"; dS = conflicting)
        stranger = "((A2.1:0.3,X9.1:0.1):0.2,(A3.1:0.2,B3.1:0.4):0.6,O2.1:1.4);"
        @test_throws ArgumentError divergence_tree(clade, "OG0000001"; dS = stranger)
        @test_throws ArgumentError divergence_tree(
            clade,
            "OG0000001";
            dS = divergence_dS,
            dN = "((A2.1:0.1,B2.1:0.1):0.1,(B3.1:0.1,A3.1:0.1):0.1,O2.1:0.1);",
        )
    end

    # -------------------------------------------------------------------------
    @testset "nearest_orthologs" begin
        codeml = mktempdir()
        write(joinpath(codeml, "OG0000001_dS.nwk"), divergence_dS)
        trees = CodemlTrees(codeml)
        by_tree = nearest_orthologs(clade, "ape")
        by_dS = nearest_orthologs(clade, "ape"; codeml = trees)
        row(table, protein) = only(eachrow(table[table.protein .== protein, :]))
        # O2.1 has both ape genes as orthologs. OrthoFinder's lengths tie (0.9 each), so
        # the first in the walk up wins; by dS, A2.1 is nearer (1.9 against 2.2).
        outgroup_tree, outgroup_dS = row(by_tree, "O2.1"), row(by_dS, "O2.1")
        @test (
            outgroup_tree.ortholog_protein,
            outgroup_tree.n_orthologs,
            outgroup_tree.basis,
        ) == ("A2.1", 2, :tree)
        @test (outgroup_dS.ortholog_protein, outgroup_dS.basis) == ("A2.1", :dS) &&
              isapprox(outgroup_dS.distance, 1.9; atol = 1e-6)
        big = row(by_dS, "B3.1")
        @test (big.gene, big.ortholog_gene, big.n_orthologs) == ("B3", "gene-SHARED-2", 1)
        # A star (OG0000003) has no lengths; a two-gene orthogroup (OG0000004) none at all.
        star = row(by_tree, "B6.1")
        @test star.basis == :none && isnan(star.distance) && star.ortholog_protein == "A5.1"
        @test !("B7.1" in by_tree.protein) && !any(by_tree.species .== "ape")
        @test ismissing(row(by_tree, "O1.1").gene)
        @test_throws ArgumentError nearest_orthologs(clade, "cat")

        # `among` leaves only the listed orthologs: O2.1 falls back to A3.1 (by gene ID,
        # dS 1.4 + 0.6 + 0.2), and B2.1, whose only ape ortholog is A2.1, drops out.
        limited = nearest_orthologs(clade, "ape"; codeml = trees, among = ["gene-SHARED-2"])
        outgroup = row(limited, "O2.1")
        @test (outgroup.ortholog_protein, outgroup.n_orthologs, outgroup.basis) ==
              ("A3.1", 1, :dS) && isapprox(outgroup.distance, 2.2; atol = 1e-6)
        @test !("B2.1" in limited.protein) && "B3.1" in limited.protein
        @test row(nearest_orthologs(clade, "ape"; among = ["A5.1"]), "B6.1").basis == :none

        # With codeml trees, ortholog pairs carry their dS: OG0000001's two ape–big
        # pairs (0.3 + 0.1 and 0.2 + 0.4), NaN for orthogroups without a tree.
        pairs = ortholog_pairs(clade, "ape", "big"; codeml = trees)
        with_dS = pairs[isfinite.(pairs.dS), :]
        @test Set(zip(with_dS.protein_a, with_dS.protein_b)) ==
              Set([("A2.1", "B2.1"), ("A3.1", "B3.1")])
        @test isapprox(sort(with_dS.dS), [0.4, 0.6]; atol = 1e-6) &&
              nrow(with_dS) < nrow(pairs)
        @test !("dS" in names(ortholog_pairs(clade, "ape", "big")))
    end

    # -------------------------------------------------------------------------
    @testset "node_ages" begin
        tree = divergence_tree(clade, "OG0000001"; dS = divergence_dS, dN = divergence_dN)
        # Mean path down to the genes below: n2 (0.3, 0.1), n3 (0.2, 0.4), n1 from
        # 0.4 deep (0.9, 0.7, 1.2, 1.4), the root over all five (1.0 for O2.1).
        @test isapprox(node_ages(tree), [1.04, 0, 0.65, 0.2, 0, 0, 0.3, 0, 0]; atol = 1e-6)
        @test node_ages(tree; rate = :dN) ≈ node_ages(tree) ./ 2
        # Pruning keeps the full tree's ages, so hiding species moves no node.
        @test isapprox(node_ages(prune(tree, "ape")), [0.65, 0, 0]; atol = 1e-6)
        # n3's deep genes outweigh n1's own mean (1.65), so n1 is raised to n3's 3.
        deep = divergence_tree(
            clade,
            "OG0000001";
            dS = "((A2.1:0.1,B2.1:0.1):0.1,(A3.1:3,B3.1:3):0.1,O2.1:0.1);",
        )
        ages = node_ages(deep)
        @test isapprox(ages[3], 3; atol = 1e-6) && ages[7] == ages[3] && ages[1] >= ages[3]
        @test all(isnan, node_ages(deep; rate = :dN))
        @test_throws ArgumentError node_ages(tree; rate = :dX)
    end
end
