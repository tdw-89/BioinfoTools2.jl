module BioinfoTools2

    const ASSETS = joinpath(pkgdir(BioinfoTools2), "assets")
    include("_load_sofa.jl")

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
