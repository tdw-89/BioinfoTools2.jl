"""
Statistical modelling over `Reference` and `Data` objects.
"""
module Modeling

using GLM
using HypothesisTests
using MultipleTesting
using StatsBase

using ..Reference
using ..Data

"""Given predictors `x`, build the design matrix of an intercept-plus-slope model."""
_design(x::AbstractVector{<:Real}) = hcat(ones(length(x)), Float64.(x))

"""
Given predictors `x` and responses `y` in `[0, 1]` — 0/1 outcomes or proportions,
such as a per-gene peak frequency — fit a logistic regression of `y` on `x`.
Returns the fitted `GLM` model; pairs with a non-finite value are dropped.

**NOTE:** a response with no variation (every `y` 0, say) has no finite fit, so it
is refused rather than fitted to a meaningless curve.
"""
function logistic_fit(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) ||
        throw(DimensionMismatch("`x` and `y` must be the same length"))
    keep = isfinite.(x) .& isfinite.(y)
    all(response -> 0 <= response <= 1, y[keep]) ||
        throw(ArgumentError("logistic responses must lie in [0, 1]"))
    allequal(y[keep]) && throw(ArgumentError("logistic responses must vary"))
    return glm(_design(x[keep]), Float64.(y[keep]), Binomial(), LogitLink())
end

"""
Given a model from [`logistic_fit`](@ref), evaluate its curve at `x`. Returns
the fitted probabilities.
"""
logistic_curve(model, x::AbstractVector{<:Real}) = predict(model, _design(x))

export logistic_fit, logistic_curve

end
