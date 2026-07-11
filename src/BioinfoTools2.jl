module BioinfoTools2

    # Internal modules
    include("so_terms.jl")

    # Exported modules
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
