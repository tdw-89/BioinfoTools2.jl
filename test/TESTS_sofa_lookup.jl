using BioinfoTools2.SOTerms
using Test

@testset "SOTermLookup" begin
    @test SO_TERMS["exon"] == SO_TERMS[:exon]
    @test SO_TERMS[SO_TERMS[147][1]][1] == 147
    @test SO_TERMS[147][1] === 0x0035
    @test SO_TERMS["gene"] == SO_TERMS[:gene]
    @test SO_TERMS[SO_TERMS[704][1]][1] == 704
    @test SO_TERMS[704][1] === 0x002e
end
