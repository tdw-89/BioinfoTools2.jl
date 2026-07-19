module BioinfoTools2

# Internal modules
include("so_terms.jl")

# Exported modules
include("reference.jl")
include("data.jl")
include("paralogs.jl")
include("plotting.jl")
include("exploration.jl")
include("modeling.jl")

using .Reference
using .Data
using .Paralogs
using .Plotting
using .Exploration
using .Modeling

export Data, Reference, Paralogs, Plotting, Exploration, Modeling
end
