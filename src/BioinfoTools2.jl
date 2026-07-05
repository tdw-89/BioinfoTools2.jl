module BioinfoTools2
    include("reference.jl")
    include("studies.jl")
    include("paralogs.jl")
    
    using .Reference
    using .Studies
    using .Paralogs

    export
        Studies,
        Reference,
        Paralogs
end
