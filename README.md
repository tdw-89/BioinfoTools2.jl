# BioinfoTools2.jl
[![codecov](https://codecov.io/gh/tdw-89/BioinfoTools2.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tdw-89/BioinfoTools2.jl)

A second attempt at creating a comprehensive suite of bioinformatics tools in pure Julia.

- [Getting Started](#getting-started)
- [Package Structure](#package-structure)
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
```

## Package Structure

| Module | What it holds |
|---|---|
| `Reference` | The in-memory genome: `Species` → `Genome` → `Scaffold` → interval tree of GFF3 features. |
| `Data` | Sample-level data loaded against a `Genome`: `BedData`, `TabularData`, `Experiment`, and interval set operations. |
| `Data.Methylation` | Single-base methylation calls, packed 8 bytes per site, with Bismark loaders and Arrow I/O. |
| `Homologs.Paralogs` | `ParalogGroup` — paralog pair relations as sparse matrices — plus reciprocal-best-hit detection. |
| `Homologs.Orthologs` | Between-genome orthology (placeholder). |
| `Exploration` | Coverage, density estimates, quantile binning and metagene profiles over the above. |
| `Modeling` | Statistical models (placeholder). |
| `Plotting` | Figures (placeholder). |

## Author
Tom Wolfe<br>
e-mail: thomas_wolfe@student.uml.edu<br>
github: [tdw-89](<https://github.com/tdw-89>)<br>