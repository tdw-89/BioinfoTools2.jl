module BioinfoTools2

# Internal modules
include("bit_codes.jl")
include("so_terms.jl")

# Exported modules
include("reference.jl")
include("data.jl")   # also defines the nested Data.Methylation submodule
include("homologs.jl")   # also defines the nested Homologs.Paralogs/Homologs.Orthologs submodules
include("plotting.jl")
include("exploration.jl")
include("modeling.jl")

using .BitCodes
using .Reference
using .Data
using .Data.Methylation
using .Homologs
using .Homologs.Paralogs
using .Homologs.Orthologs
using .Plotting
using .Exploration
using .Modeling

# `Methylation` is nested inside `Data`, and `Paralogs`/`Orthologs` inside
# `Homologs`, but all three are re-exported here as well so they can still be
# reached (and `using`d) as `BioinfoTools2.Methylation` etc.
export BitCodes,
    Data,
    Reference,
    Methylation,
    Homologs,
    Paralogs,
    Orthologs,
    Plotting,
    Exploration,
    Modeling

end
