using BioinfoTools2.Modeling
using Test

@testset "Modeling" begin
    @testset "logistic_fit and logistic_curve" begin
        # Proportions lying exactly on a known curve, logit(p) = -2 + 0.004x,
        # recover its coefficients, as a peak frequency along a ranking would.
        x = collect(0.0:50.0:2000.0)
        p = @. 1 / (1 + exp(-(-2 + 0.004 * x)))
        model = logistic_fit(x, p)
        @test Modeling.coef(model) ≈ [-2.0, 0.004] rtol = 1e-4
        @test logistic_curve(model, [0.0, 500.0]) ≈ p[[1, 11]] rtol = 1e-4

        # Non-finite pairs are dropped rather than failing the fit.
        with_gap = logistic_fit([x; NaN], [p; 0.5])
        @test Modeling.coef(with_gap) ≈ Modeling.coef(model)

        @test_throws ArgumentError logistic_fit([1.0, 2.0], [0.5, 1.5])
        @test_throws DimensionMismatch logistic_fit([1.0], [0.5, 0.5])
        # All-zero responses have no finite fit (the slope runs off to -Inf).
        @test_throws ArgumentError logistic_fit([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
    end
end
