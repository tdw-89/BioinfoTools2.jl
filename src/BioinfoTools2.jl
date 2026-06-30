module BioinfoTools2
    include("reference.jl")
    include("studies.jl")
    include("paralogs.jl")
    using .Paralogs

    export Paralogs
end
