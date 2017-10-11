mutable struct ScaledFit{T, fitType}
    μ::T
    σ::T
    fit_func::fitType
end
ScaledFit(μ, σ, fit_func) = ScaledFit{typeof(μ), typeof(fit_func)}(
                                μ, σ, fit_func)

function (f::ScaledFit)(x) f.fit_func((x-f.μ)/f.σ) end
function der(f::ScaledFit, query, order=1)
    der(f.fit_func, order=order)((query-f.μ)/f.σ)/f.σ^order
end
"""
Alias for polyder when fit function is polynomial
"""
function der(f::Poly, order=1) polyder(f, order) end
