"""
Homology relations between genes.

A thin parent over two submodules, split by the kind of homology rather than by
the machinery: [`Paralogs`](@ref) for duplicates *within* a genome and
[`Orthologs`](@ref) for correspondences *between* genomes. Both are re-exported
here, so `BioinfoTools2.Paralogs` reaches the same module as
`BioinfoTools2.Homologs.Paralogs`.
"""
module Homologs

include("homologs/paralogs.jl")
include("homologs/orthologs.jl")

using .Paralogs
using .Orthologs

export Paralogs, Orthologs

end
