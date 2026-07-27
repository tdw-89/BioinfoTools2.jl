module BioinfoTools2

# Internal modules
include("bit_codes.jl")
include("so_terms.jl")

# Exported modules
include("reference.jl")
include("data.jl")
include("methylation.jl")
include("paralogs.jl")
include("plotting.jl")
include("exploration.jl")
include("modeling.jl")

using .BitCodes
using .Reference
using .Data
using .Methylation
using .Paralogs
using .Plotting
using .Exploration
using .Modeling

export BitCodes, Data, Reference, Methylation, Paralogs, Plotting, Exploration, Modeling
end
