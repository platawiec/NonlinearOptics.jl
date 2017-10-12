mutable struct ScaledFit{T, fitType}
    μ::T
    σ::T
    fit_func::fitType
end
ScaledFit(μ, σ, fit_func) = ScaledFit{typeof(μ), typeof(fit_func)}(
                                μ, σ, fit_func)

function (f::ScaledFit)(x) f.fit_func((x-f.μ)/f.σ) end
function der(f::ScaledFit, query; order=1)
    der(f.fit_func, order=order)((query-f.μ)/f.σ)/f.σ^order
end
"""
Alias for polyder when fit function is polynomial
"""
function der(f::Poly; order=1) polyder(f, order) end
function der(f::Poly, query; order=1) polyder(f, order)(query) end
"""
Macro for call work-around See JuliaLang Pull #23168
https://github.com/JuliaLang/julia/pull/23168
"""
t_info(ex::Symbol) = (ex, tuple())
t_info(ex::Expr) = ex.head == :(<:) ? t_info(ex.args[1]) : (ex, ex.args[2:end])
macro ev(fsig)
    @assert fsig.head == :type
    tname, tparams = t_info(fsig.args[2])
    sym = esc(:x)
    return quote
        $(esc(fsig))
        ($sym::$tname)(args...) = prop_call($sym, args...)
    end
end
