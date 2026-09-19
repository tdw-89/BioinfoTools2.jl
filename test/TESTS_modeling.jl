using BioinfoTools2.Modeling
using Test

@testset "Modeling" begin
    @testset "binarize" begin
        labels = binarize(Dict("a" => 0.9, "b" => 0.5, "c" => 0.1, "d" => NaN), 0.5)
        @test labels == Dict("a" => 1.0, "b" => 1.0, "c" => 0.0)
        # An unmeasured gene is dropped, not called negative.
        @test !haskey(labels, "d")
    end

    @testset "logistic_fit and logistic_curve" begin
        # Labels drawn deterministically from logit(p) = -2 + 0.004x: every
        # predictor left of the boundary is negative, every one right positive.
        x = collect(0.0:50.0:2000.0)
        y = Float64.(x .>= 500.0)
        model = logistic_fit(x, y)
        curve = logistic_curve(model, x)
        @test all(0 .<= curve .<= 1)
        @test issorted(curve)
        @test curve[1] < 0.5 < curve[end]

        # Non-finite pairs are dropped rather than failing the fit.
        with_gap = logistic_fit([x; NaN], [y; 1.0])
        @test Modeling.coef(with_gap) ≈ Modeling.coef(model)

        # Proportions are no longer a valid response: this is a classifier.
        @test_throws ArgumentError logistic_fit([1.0, 2.0], [0.25, 0.75])
        @test_throws ArgumentError logistic_fit([1.0, 2.0], [0.0, 1.5])
        @test_throws DimensionMismatch logistic_fit([1.0], [0.0, 1.0])
        # One class alone has no finite fit (the slope runs off to -Inf).
        @test_throws ArgumentError logistic_fit([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
    end

    @testset "decision_boundary" begin
        # A clean split at 10 puts the crossing between the two classes.
        x = collect(1.0:20.0)
        model = logistic_fit(x, Float64.(x .> 10))
        boundary = decision_boundary(model)
        @test 10 < boundary < 11
        @test logistic_curve(model, [boundary]) ≈ [0.5]

        # Labels independent of the predictor never cross, so the boundary is
        # pushed off the axis rather than landing somewhere arbitrary.
        flat = decision_boundary(logistic_fit([1.0, 1.0, 2.0, 2.0], [0.0, 1.0, 0.0, 1.0]))
        @test !isfinite(flat) || abs(flat) > 1e6
    end
end
