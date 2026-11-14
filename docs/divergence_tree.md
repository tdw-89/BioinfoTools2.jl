# `DivergenceTree` reference

How `Homologs.Orthologs.DivergenceTree` puts codeml's per-branch dS and dN onto
an orthogroup of a [`Clade`](clade.md), the rules it follows, how
`Plotting.gene_family_plot!` draws it, and how it checks out on the two
divergence-pipeline runs.

## Building one

```julia
using BioinfoTools2.Orthologs

trees = CodemlTrees("results/codeml/trees")      # <OG>_dS.nwk (+ <OG>_dN.nwk)
tree = divergence_tree(clade, trees, "OG0000100")
tree = divergence_tree(clade, "OG0000100"; dS = newick, dN = nothing)   # from strings
```

`CodemlTrees(dir)` lists the orthogroups with a `<OG>_dS.nwk`, sorted; a
`<OG>_dN.nwk` beside it is read too (it must have the same shape and labels),
else `dN` is `NaN`. The divergence pipeline writes one per orthogroup with
paralogs: 8,945 of the fish run's 17,740 resolved trees, 2,842 of the primates'.

## Input

codeml's trees are the gene tree it was fitted on, with each branch's dS (or
dN) as its length:

- **Leaves** are bare protein accessions (`XP_069040797.1`), matched to the
  orthogroup's forest leaves by `GeneForest.protein`. An accession the
  orthogroup lacks, or one leaf twice, throws.
- **Unrooted**: written as a basal trifurcation (or, for a few, a bifurcation,
  whose two branches are one).
- **Some genes are missing** (`dropped`): codeml ran on a filtered alignment.
- **Polytomies of OrthoFinder's tree are resolved**, sometimes with long
  branches (up to dS 6 in the fish run). These **codeml-only splits** have
  `exact == false`.

The pairwise `tree_dS` / `tree_dN` columns of `<OG>_ks.tsv` are these trees'
path sums, so `patristic(tree, a, b)` reproduces them.

## Rooting

codeml's tree is drawn in its own topology, rooted where OrthoFinder roots the
orthogroup:

1. Take the forest LCA of the leaves codeml kept, and group those leaves by the
   LCA's children (OrthoFinder's root clades, as far as codeml has them).
2. **Two groups:** they are split by one codeml branch. The root goes on it, the
   branch's dS (and dN) divided in proportion to OrthoFinder's own lengths from
   each group's LCA up to the root (half each when those are missing or zero).
3. **Three or more** (a root polytomy): the codeml nodes outside every group are
   where OrthoFinder's polytomy falls. The root is the one whose farthest group
   is nearest by dS (ties to the first in codeml's order).
4. Children are ordered by their first forest leaf, so leaves follow
   OrthoFinder's Newick order.

The root's `dS` is `NaN`. Writing the same unrooted tree rooted elsewhere gives
the same `DivergenceTree`.

## Mapping onto the forest

- `forest_node[node]` is a leaf's forest leaf, else the forest LCA of its leaves.
- `exact[node]`: the forest node holds exactly this node's codeml leaves, i.e.
  OrthoFinder has this clade. A codeml-only split shares its parent's
  `forest_node`.
- **Compatibility** is checked: a codeml-only clade must be a union of whole
  children of its forest node (counting kept leaves). A tree that splits an
  OrthoFinder clade throws.
- `species[node]` is the node's **lineage** mask: the species under the child of
  its parent's forest node that holds it, so genes codeml dropped still count.
  At a codeml-only split it is the OR of its children's; at the root, the forest
  node's own.

## Events

`event[node]` applies the per-child-pair rule of [`clade.md`](clade.md) to the
children's lineage masks: `:duplication` when every pair shares a species,
`:speciation` when none does, `:mixed` otherwise; `:leaf` at a leaf, and
`:unresolved` throughout a `TREE_PAIR` orthogroup. On a node whose children are
all OrthoFinder's, this is OrthoFinder's call (`node_event(::GeneTree, …)`) even
when codeml dropped genes. At a codeml-only split it reconciles codeml's
resolution by species overlap.

`species_node(tree, node)` places a node on the species tree as OrthoFinder
does: by its OrthoFinder clade's species (its lineage mask at a codeml-only
split); a leaf's is its species.

## Pruning

`prune(tree, species)` keeps the named species' leaves:
- A node left with one child merges into one branch, lengths summed.
- When the root goes, the new root's branch is the **stem** down from the full
  tree's root, so pruned trees share its time axis. `root_distances` starts at
  the new root. The plot draws the stem to the left of 0, but a stem longer than
  a tenth of the tree's depth runs off the axis edge instead of squeezing the
  tree; its full length is in the root's tooltip.
- Each kept node takes the lineage mask at the top of its merged chain, and
  events are re-decided over the surviving children with those masks. A
  duplication whose other copy was pruned therefore stays a duplication.
- Leaves are kept by their own species, not their lineage mask.
- `dropped` is filtered to the kept species.

## Ages

`node_ages(tree; rate = :dS)` dates every node by **mean path length**: the mean
`rate` from the node down to all the leaves below it. Leaves are 0. A parent younger
than one of its children is raised to that child's age, so ages never decrease
towards the root. For a duplication of two single genes, the age is half their
pairwise dS (the usual "Ks age").

Ages are computed once, on the unpruned tree, and `prune` keeps them, so hiding
species moves no node. They assume a clock, so a fast lineage pulls up the age
of every node above it. In primate OG0000363, one macaque gene's dS 4.1 branch
dates the root at ~2.1, while the rest of the tree is under 0.2. The clamp then
passes such an age further up. `dN` ages are `NaN` without a dN tree.

## Nearest orthologs

`nearest_orthologs(clade, species; codeml = trees)` gives each gene of the other
species its nearest ortholog in `species`, among the orthologs the per-child-pair
rule calls (so none in a `TREE_PAIR` orthogroup):
- **`:dS`** when the orthogroup's divergence tree holds the gene and every one of
  those orthologs: the smallest `patristic` dS.
- **`:tree`** otherwise (no codeml tree, or codeml dropped one of them): the
  smallest path length on OrthoFinder's tree. Its lengths are protein
  substitutions per site, so `distance` is only comparable within a `basis`.
- **`:none`** in a star (`TREE_STAR`), which has no lengths: the first ortholog.

Ties go to the ortholog met first walking up from the gene. `n_orthologs` counts
the candidates; a gene with one is unaffected by the metric. Rows are keyed by
`gene` (`missing` for a genome-less species) and `protein`.

`among` (gene IDs or protein accessions of `species`) drops every other candidate
before the nearest is chosen, so each gene gets its nearest ortholog *within* that
set (e.g. the nearest one with a score), `n_orthologs` counts only those, and a gene
with none there gets no row. The basis is then decided over the remaining
candidates alone.

`ortholog_pairs(clade, a, b; codeml = trees)` adds a `dS` column to the pair table: each pair's `patristic` dS on its orthogroup's divergence tree, `NaN` without a codeml tree or where it lacks either gene.

## Drawing

`gene_family_plot!(position, tree; layout, tracks, rate = :dS, saturation = 2.0, …)`
draws the tree with one of two `layout`s:
- **`:chronogram`** (the default): every gene at age 0 on the right and each
  node at its `node_ages`, older to the left. Each branch runs from its parent's
  age to its own, coloured by its `relative_rate`: log2 of its `rate` over the age
  it spans, from ¼× (blue) through 1× (grey) to 4× (red), clamped. This keeps
  visible the rate differences the alignment hides. A pruned tree's stem is a
  short stub.
- **`:phylogram`**: x is `rate` summed from the root, so the distance between two
  genes along the tree is their pairwise dS.

In both, branches longer than `saturation` (dS saturates around 2–3) are dashed.
The returned `positions` are the nodes' `(x, y)`.
- **Leaves**: coloured by `species_palette(clade)`, one fixed colour per species
  of the clade.
- **Internal nodes**, by event: red diamonds for duplications, open circles for
  speciations, orange squares for mixed polytomies. Codeml-only splits are faded.
- **Tracks** are heatmap columns beside the leaves:
  - a `GeneTrack(name, values; colormap, colorrange, label, notes)`. `values` is
    `species name => (gene ID or protein accession => value)`, per species what
    `Exploration.tss_window` returns. A leaf is looked up by gene ID, then by
    protein, which is how genome-less species get values;
  - or a plain `name => values` pair, drawn with the plot's `colormap`,
    `colorrange` and `track_label`.

  Each track keeps its own colour scale. Tracks sharing `label`, colormap and
  range share one colorbar; the colorbars sit side by side, followed by the
  branch-rate scale in a chronogram. Missing values are grey,
  and `notes` (the same shape, holding strings) are added to the tooltips. `opacity` (the same shape, in [0, 1]) fades each cell over the white background; a gene without one is opaque.
  `why_missing` (the same shape, holding strings) gives a missing cell a reason,
  drawn in its colour from `missing_colors` (`reason => color` pairs; a reason
  without one stays grey) and shown in the cell's tooltip. Each track with
  `missing_colors` adds them to a "Missing track values" legend group, after the
  grey, which is labelled `missing_label` (`"no value"` by default).
- **Labels** are gene IDs without `gene-`, else proteins, as right-hand ticks.
  They are hidden past `max_labels` leaves.
- **Tooltips** (`inspector_label`, shown under a `DataInspector`):
  - leaves: gene, protein, species and track values;
  - internal nodes: event, species-tree node, whether OrthoFinder or codeml made
    the split, genes below, and dS age;
  - every node: branch dS and dN, and distance from the root.

`gene_family_explorer(clade, trees; species, tracks, …)` is a GLMakie window
over the same plot. It lives in `ext/BioinfoTools2GLMakieExt.jl`, which loads on
`using GLMakie`.
- A textbox takes an orthogroup, gene ID or protein; a searched gene's leaf is
  ringed.
- A checkbox per species of the clade prunes every tree to the checked species
  (`species` sets the initial ones; all by default). Changing them re-selects
  the orthogroups and stays on the current one while it still qualifies.
- An "align genes (ages)" toggle switches between the chronogram (on, by
  default; `layout`) and the phylogram, redrawing the current tree.
- ◀/▶ step through the orthogroups that have a codeml tree and two or more genes
  of the checked species.
- Hovering shows the tooltips.

## Agreement (OrthoFinder 3.1.3 + codeml, both divergence runs)

Checked by building every tree, with no genomes loaded:

| Check | fish (8 species) | primates (7 species) |
|---|---|---|
| trees built (none conflicts with OrthoFinder) | 8,945 of 8,945 | 2,842 of 2,842 |
| `patristic` vs `<OG>_ks.tsv` `tree_dS`, max error | 612,657 pairs, 3.3e-6 | 130,272 pairs, 1.0e-4 (4-decimal rounding) |
| same for `tree_dN` | 1.1e-7 | 2.0e-7 |
| events = `node_event(::GeneTree)` where children are OrthoFinder's | 102,987 of 102,987 | 29,800 of 29,800 |
| codeml-only splits | 1,129 | 3,021 |
| genes codeml dropped | 862 | 1,419 |

Building every fish tree takes ~2.3 s (0.26 ms each; 0.7 ms for the
690-gene OG0000000).
