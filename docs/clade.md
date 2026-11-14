# `Clade` reference

How `Homologs.Orthologs.Clade` stores a group of species and their OrthoFinder
gene trees, the rules its queries follow, and how closely they reproduce
OrthoFinder's own tables.

## Building one

```julia
using BioinfoTools2.Orthologs

clade = load_clade(
    "results/orthofinder",                       # an OrthoFinder results directory
    Dict(
        "human" => "human/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff",
        "chimp" => "chimp/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri_genomic.gff",
        # ...one GFF3 per species you have a genome for
    ),
)
```

`load_clade` loads each GFF3 with `add_features!` (genes only by default, set by
`load_options`) and reads its `protein_gene_map`, one task per species. It then
calls `Clade(species, dir; protein_maps)`, which also accepts already-loaded
`Species`. The keys are OrthoFinder's species names, which are its input file
names and the leaves of its species tree.

A species of the tree with no GFF3 (in the primate run, mouse and bonobo) gets
an **empty genome**. Its genes stay in the trees and still decide every call,
which the outgroup has to. It is queried by protein accession, and its `gene`
columns are `missing`.

## Input files

Read from the results directory, each path overridable by keyword:

| Keyword | Default | Holds |
|---|---|---|
| `species_tree` | `Species_Tree/SpeciesTree_rooted_node_labels.txt` | the rooted species tree |
| `gene_trees` | `Resolved_Gene_Trees/Resolved_Gene_Trees.txt` | `OG…: <newick>` per orthogroup with a tree |
| `orthogroups` | `Orthogroups/Orthogroups.tsv` | every orthogroup's members (`nothing` skips) |

What the readers expect, and throw on otherwise:

- **Gene-tree leaves** are `<species>_<protein>`. They are split on the
  *longest* species name that prefixes them, since a name may hold `_`
  (`pongo_abelii_XP_054409609.1`).
- **Internal gene-tree labels** are OrthoFinder's `n0, n1, …` in preorder (the
  root is `n0`). They aren't stored; `node_label` recomputes them. Empty or
  support-valued labels are tolerated.
- **Species-tree internal nodes** without a name, or named by a support value
  (`SpeciesTree_rooted.txt`), are named `N0, N1, …` in preorder, the same names
  `SpeciesTree_rooted_node_labels.txt` gives.
- At most 64 species.

## Layout

All trees live in one `GeneForest`, concatenated in **preorder**, with nodes
indexed by `Int32`. Struct-of-arrays, no node objects:

| Field | Per | Meaning |
|---|---|---|
| `parent::Vector{Int32}` | node | parent node, 0 at a root |
| `subtree_end::Vector{Int32}` | node | last descendant: the subtree is `node:subtree_end[node]` |
| `species::Vector{UInt64}` | node | bitmask of the species below; bit `s - 1` is the `s`-th species |
| `branch_length::Vector{Float32}` | node | `NaN` where absent (roots, star trees) |
| `gene::Vector{UInt32}` | node | a leaf's metadata index in its species' genome; 0 if none |
| `protein::Vector{String}` | node | a leaf's protein accession; `""` for internal nodes |
| `tree_start::Vector{Int32}` | tree | first node of each tree (`n_trees + 1` entries) |
| `orthogroups::Vector{String}` | tree | names, **sorted**, so a name is found by binary search |
| `source::Vector{UInt8}` | tree | `TREE_RESOLVED`, `TREE_STAR` or `TREE_PAIR` |

What preorder buys, with no child lists:
- **Leaf test:** a node is a leaf iff `subtree_end[node] == node`.
- **Children:** `node + 1`, then each `subtree_end[child] + 1`, up to `subtree_end[node]`.
- **Descendant test:** "is `y` below `a`" is `a <= y <= subtree_end[a]`.
- **Leaves under a node:** one contiguous range.

`SpeciesTree` uses the same layout for its one tree. It adds `labels`, a
`depth` (root 1, as `ParalogGroup.lca_depth` counts) and a mask per node.

Each species also has two `GeneIndex`es, by gene ID and by protein accession.
These are sorted keys searched by binary search, not a `Dict`, the way
`Reference` resolves IDs. Gene IDs are the genome's own vocab strings, so the
index adds no copies.

Measured on the primate run (7 species, 19,927 orthogroups, 270,532 nodes):
- forest: 12.7 MiB;
- gene and protein indexes: 8.0 MiB;
- the five genes-only genomes: ~4.5 MiB each.

## Rules

**Orthologs.** Two leaves are orthologs iff the two children of their common
ancestor that hold them share no species. Implementations:
- `orthologs` walks up from the leaf, taking every sibling disjoint from the
  side it came from.
- `ortholog_pairs` checks each pair of children at each node.

Deciding per *child pair* is what handles polytomies. A polytomy can hold pairs
of both kinds, and then it is neither a duplication nor a speciation as a whole.
OrthoFinder's `Duplications.tsv` lists a node only when *every* pair of its
children overlaps: all such binary nodes, plus polytomies that overlap
throughout. It leaves out the mixed polytomies.

**Events.** `common_ancestor` reports `:duplication` when those two children
share a species, and `:speciation` otherwise.

**Placement.** A gene-tree node maps to the deepest species-tree node holding
all its species. `depth` is that node's.

**Tree-less orthogroups.** `Orthogroups.tsv` orthogroups with no resolved tree
come in as stars, one root over every member:
- **Three or more genes** (`TREE_STAR`): every cross-species pair is called.
- **Two genes** (`TREE_PAIR`): no orthologs are called, and `common_ancestor`
  reports `:unresolved`, because OrthoFinder 3.1.3 calls none there. These
  pairs are probably real one-to-one orthologs, but OrthoFinder's tables are
  kept as the reference.

**`one_to_one`.** A pair is one-to-one when each gene is the other's only
ortholog in the other's species.

**Genes.** `protein_gene_map` follows each CDS's `protein_id` up its `Parent`
chain (CDS → mRNA → gene, or CDS → gene) to a `gene` record. This is exact
where `"gene-" * symbol` is not (`gene-TRNAK-UUU-2`). A leaf whose protein maps
to no gene, or to a gene not loaded, keeps gene 0, and the build warns with
both counts.

## Agreement with OrthoFinder (primate run, OrthoFinder 3.1.3)

| Check | Result |
|---|---|
| `ortholog_pairs` vs `Orthologues/<a>__v__<b>.tsv`, every species pair | 21 of 21 equal (525,277 pairs) |
| `Duplications.tsv` rows: event, species-tree node and gene-tree node of the first `Genes 1 × Genes 2` pair | 12,196 of 12,196 (266 rows list no genes on one side) |
| `Duplications.tsv` nodes: every child pair overlapping, placed at the listed node | 12,462 of 12,462 (12,196 binary, 266 polytomies) |
| `newick` round trip of `Resolved_Gene_Trees.txt` (lengths at `Float32`) | 18,701 of 18,701 |
| leaves placed on a gene, the five species with a GFF3 | 106,095 of 106,095 |

## Costs

| Operation | Time (primate run, 8 threads) |
|---|---|
| `load_clade`, five NCBI GFF3s | ~7 s (each GFF3 is read twice: genes, then protein links) |
| `orthologs`, per gene | ~3.4 µs |
| `orthogroup`, per gene | ~0.4 µs |
| `ortholog_pairs`, all 21 species pairs | ~0.9 s |

## codeml branch lengths

`divergence_tree` roots a divergence-pipeline codeml tree (per-branch dS and dN)
on an orthogroup of the clade, for `Plotting.gene_family_plot!`. See
[`divergence_tree.md`](divergence_tree.md).

## Not carried

- **HOGs** (`Phylogenetic_Hierarchical_Orthogroups`); they can be derived from
  the forest and its masks.
- **`Duplications.tsv` support values**; events come from the trees.
- **Unassigned genes** (`Orthogroups_UnassignedGenes.tsv`), which have no
  orthogroup.
- **Several leaves for one gene** (several isoforms in the proteome) build, but
  warn, and the first leaf answers lookups.
