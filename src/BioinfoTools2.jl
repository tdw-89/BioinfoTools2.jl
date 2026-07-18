module BioinfoTools2

# Internal modules
include("so_terms.jl")

# Exported modules
include("reference.jl")
include("studies.jl")
include("paralogs.jl")
include("plotting.jl")
include("exploration.jl")
include("modeling.jl")

using .Reference
using .Studies
using .Paralogs
using .Plotting
using .Exploration
using .Modeling

export Studies, Reference, Paralogs, Plotting, Exploration, Modeling
end
