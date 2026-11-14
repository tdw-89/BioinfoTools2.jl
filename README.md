# BioinfoTools2.jl
[![codecov](https://codecov.io/gh/tdw-89/BioinfoTools2.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tdw-89/BioinfoTools2.jl)

A second attempt at creating a comprehensive suite of bioinformatics tools in pure Julia.

- [Getting Started](#getting-started)
- [Package Structure](#package-structure)
- [Reference documents](#reference-documents)
- [Author](#author)

## Getting Started

The package depends on forked versions of several BioJulia packages that are not
in the public registry, so add those before installing it:

```julia
using Pkg
Pkg.add([
    PackageSpec(url="https://github.com/tdw-89/Indexes.jl"),
    PackageSpec(url="https://github.com/tdw-89/GenomicFeatures.jl"),
    PackageSpec(url="https://github.com/tdw-89/GFF3.jl.git"),
    PackageSpec(url="https://github.com/tdw-89/BED.jl.git"),
])
Pkg.add(url="https://github.com/tdw-89/BioinfoTools2.jl")
```

To run the tests from a clone, unpack the fixtures first — only the tarball is
tracked:

```sh
tar -xzf test/data/test_data.tar.gz -C test/data
julia --project=. -t auto -e 'using Pkg; Pkg.test()'
```

## Package Structure

| Module | What it holds |
|---|---|
| `BitCodes` | The bit-packing vocabulary the other modules share: strand codes and packed-field helpers. |
| `Reference` | The in-memory genome: `Species` → `Genome` → `Scaffold` → interval tree of GFF3 features. |
| `Data` | Sample-level data loaded against a `Genome`: `BedData`, `TabularData`, `Experiment`, and interval set operations. |
| `Data.Methylation` | Single-base methylation calls, packed 8 bytes per site, with Bismark loaders and Arrow I/O. |
| `Data.Mappability` | Per-base empirical mappability from [UniMap.jl](https://github.com/tdw-89/UniMap.jl): converts a `mappability.zarr` (with `using Zarr`) into a lossless run-length Arrow cache, read lazily per channel, and used for metagenes, MAPQ-weighted recall, low-mappability intervals and masking methylation calls. |
| `Homologs.Paralogs` | `ParalogGroup` — paralog pair relations as sparse matrices — plus divergence-table loading, reciprocal-best-hit detection, duplication-node annotation from OrthoFinder output, and ranking genes by duplication age. |
| `Homologs.Orthologs` | `Clade`: several `Species`, their species tree, and every OrthoFinder gene tree in one compact preorder forest, for ortholog lookups, pairwise ortholog tables, common ancestors (duplication vs speciation, and where in the species tree) and gene-tree walks; `DivergenceTree` roots codeml's per-branch dS/dN trees on those gene trees, and `nearest_orthologs` finds each gene's closest ortholog in one species. |
| `Exploration` | Coverage, density estimates, quantile binning, metagene profiles, per-gene profile matrices and TSS-window summaries over the above. |
| `Modeling` | Statistical models: logistic regression as a binary classifier. |
| `Plotting` | Metagene line profiles and heatmaps, rank heatmaps, per-quantile violin and box plots, rank-trend plots, fitted probability curves, and gene-family phylograms with dS as distance and per-gene heatmap tracks (CairoMakie; an interactive explorer with `using GLMakie`). |

## Reference documents

- [`docs/paralog_group.md`](docs/paralog_group.md) — `ParalogGroup`'s input columns, relation-matrix layouts, and `lca_depth` numbering.
- [`docs/weighting.md`](docs/weighting.md) — how read depth becomes weight through the methylation metagene pipeline.
- [`docs/mappability.md`](docs/mappability.md) — converting UniMap output, the cache layout and sizes, genome matching, and the analysis built on it.
- [`docs/clade.md`](docs/clade.md) — building a `Clade` from an OrthoFinder run, the gene-forest layout, the ortholog and placement rules, and how they reproduce OrthoFinder's tables.
- [`docs/divergence_tree.md`](docs/divergence_tree.md) — rooting codeml's dS/dN trees on a `Clade`'s gene trees, the event rules, pruning, nearest orthologs, the gene-family plot and explorer, and the checks against both divergence runs.

## Author
Tom Wolfe<br>
e-mail: thomas_wolfe@student.uml.edu<br>
github: [tdw-89](<https://github.com/tdw-89>)<br>