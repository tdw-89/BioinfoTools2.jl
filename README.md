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
| `Homologs.Paralogs` | `ParalogGroup` — paralog pair relations as sparse matrices — plus reciprocal-best-hit detection and duplication-node annotation from OrthoFinder output. |
| `Homologs.Orthologs` | Between-genome orthology (placeholder). |
| `Exploration` | Coverage, density estimates, quantile binning, metagene profiles and TSS-window summaries over the above. |
| `Modeling` | Statistical models: logistic regression as a binary classifier. |
| `Plotting` | Metagene line profiles and heatmaps, per-quantile violin plots, rank-trend plots, and fitted probability curves (CairoMakie). |

## Reference documents

- [`docs/paralog_group.md`](docs/paralog_group.md) — `ParalogGroup`'s input columns, relation-matrix layouts, and `lca_depth` numbering.
- [`docs/weighting.md`](docs/weighting.md) — how read depth becomes weight through the methylation metagene pipeline.

## Author
Tom Wolfe<br>
e-mail: thomas_wolfe@student.uml.edu<br>
github: [tdw-89](<https://github.com/tdw-89>)<br>