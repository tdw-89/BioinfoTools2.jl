using BioinfoTools2.Reference

const DATA_DIR = joinpath(@__DIR__, "data")
const GFF_SINGLE = joinpath(DATA_DIR, "NC_003280.10.gff.gz")
const REFERENCE_SUITE = BenchmarkGroup()

## add_features!
REFERENCE_SUITE["add_features! - Single"] =
    @benchmarkable add_features!(gff_file_path, test_genome) setup(
        gff_file_path = GFF_SINGLE,
        test_genome = Species("C. elegans").genome,
    )
