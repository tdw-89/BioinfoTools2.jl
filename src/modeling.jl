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
Given a `key => value` `Dict` and a `threshold`, label each value `1.0` at or
above it and `0.0` below. Returns a `Dict` over the same keys, ready to be
ranked and handed to [`logistic_fit`](@ref).

**NOTE:** a non-finite value is unmeasured, not negative, so its key is dropped
rather than labelled `0`.
"""
binarize(values::AbstractDict{K,<:Real}, threshold::Real) where {K} = Dict{K,Float64}(
    key => Float64(value >= threshold) for (key, value) in values if isfinite(value)
)

"""
Given predictors `x` and class labels `y`, each `0` or `1`, fit a logistic
regression of P(`y` = 1) on `x`. Returns the fitted `GLM` model; pairs with a
non-finite value are dropped.

**NOTE:** labels of one class only have no finite fit — the slope runs off to
±Inf — so they are refused rather than fitted to a meaningless curve.
"""
function logistic_fit(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) ||
        throw(DimensionMismatch("`x` and `y` must be the same length"))
    keep = isfinite.(x) .& isfinite.(y)
    labels = Float64.(y[keep])
    all(label -> label == 0 || label == 1, labels) ||
        throw(ArgumentError("logistic labels must each be 0 or 1"))
    allequal(labels) && throw(ArgumentError("logistic labels must cover both classes"))
    return glm(_design(x[keep]), labels, Binomial(), LogitLink())
end

"""
Given a model from [`logistic_fit`](@ref), evaluate its curve at `x`. Returns
the predicted probability of the positive class at each point.
"""
logistic_curve(model, x::AbstractVector{<:Real}) = predict(model, _design(x))

"""
Given a model from [`logistic_fit`](@ref), locate the predictor value at which
it predicts either class equally. Returns `-intercept / slope`, the `x` where
the curve crosses `0.5`, which is non-finite when the slope is flat.

**NOTE:** the model is unbounded, so the boundary may fall outside the fitted
range of `x`, meaning the classifier never switches over the data — check
`isfinite` and the range before drawing it.
"""
function decision_boundary(model)
    intercept, slope = coef(model)
    return -intercept / slope
end

export binarize, logistic_fit, logistic_curve, decision_boundary

end
