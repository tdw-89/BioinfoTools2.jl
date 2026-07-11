using BioinfoTools2
using Test

@testset "BioinfoTools.jl" begin
    include("TESTS_sofa_lookup.jl")
    include("TESTS_paralogs.jl")
    include("TESTS_reference.jl")
    include("TESTS_studies.jl")
end