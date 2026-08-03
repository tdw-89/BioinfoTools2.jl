module BioinfoTools2

# Internal modules
include("bit_codes.jl")
include("so_terms.jl")

# Exported modules
include("reference.jl")
include("data.jl")   # also defines the nested Data.Methylation submodule
include("paralogs.jl")
include("plotting.jl")
include("exploration.jl")
include("modeling.jl")

using .BitCodes
using .Reference
using .Data
using .Data.Methylation
using .Paralogs
using .Plotting
using .Exploration
using .Modeling

# `Methylation` is nested inside `Data`, but re-exported here as well so it can
# still be reached (and `using`d) as `BioinfoTools2.Methylation`.
export BitCodes, Data, Reference, Methylation, Paralogs, Plotting, Exploration, Modeling
end
