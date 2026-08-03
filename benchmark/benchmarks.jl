using BenchmarkTools
using BioinfoTools2

# Top-level benchmark suite for the package
const SUITE = BenchmarkGroup()

# Public modules

## Reference
include("BENCH_Reference.jl")
SUITE["Reference"] = REFERENCE_SUITE

## Data
include("BENCH_Data.jl")

## Data.Methylation
include("BENCH_Data_Methylation.jl")

## Paralogs

## Exploration

# Internal modules

## SOTerms

## BitCodes