# BioinfoTools2.jl
[![codecov](https://codecov.io/gh/tdw-89/BioinfoTools2.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tdw-89/BioinfoTools2.jl)

A second attempt at creating a comprehensive suite of bioinformatics tools in pure Julia.

- [Getting Started](#getting-started)
- [Package Structure](#package-structure)
- [Author](#author)

## Getting Started

## Package Structure
As the purpose of this package is to allow for the easy organization, analysis and manipulation of bioinformatics data, every type unique to this package is organized hierarchically in a [*Study*](<./src/studies.jl>):
```mermaid
classDiagram
direction TB
    class Study {
	    String id
	    String title
	    Date date
    }

    class Assay {
	    String id
	    String type
	    String description
    }

    class Measurement {
	    String file_path
	    String format
	    Float data_size
    }

    class Data {
        <<Union>>
    }

    class Tabular {
        Vector~String~ variables
        Vector~Tuple~String, UInt32~~ samples
        Matrix table
    }

    class BedData {
        Dict~String, IntervalMeta64~  scaffolds
    }

    Data <|-- Tabular
    Data <|-- BedData

    class BioSample {
	    String sample_id
	    String tissue_type
	    Species species
    }

    class Species {
	    String name
	    String taxon_id
    }

    class Genome {
        Dict~String, Scaffold~ scaffolds
    }

    class AssayMethod {
	    String name
	    String description
    }

    Study "1" --> "*" Assay
    Assay --> Measurement
    Measurement --> Data
    Assay --> BioSample
    BioSample --> Species
    Species --> Genome
    Data .. Genome
    Assay --> AssayMethod
    AssayMethod .. Measurement
```
As shown, a `Study` is composed of one or more `Assay`'s. An `Assay` contains a `Measurement` with processed data (currently either BED/interval-based data or tabular data).


## Author
Tom Wolfe<br>
e-mail: thomas_wolfe@student.uml.edu<br>
github: [tdw-89](<https://github.com/tdw-89>)<br>