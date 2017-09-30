"""
∂A/∂t = N A + D A
D is dispersive operator, N is nonlinear operator
Note that for NLSE for pulse propagation, substitute z -> t
"""
mutable struct NLSEProblem{uType, tType, zType, F1, F2, C} <: AbstractNLSEProblem{uType, tType, zType}
    N::F1
    D::F2
    u0::uType
    tspan::Tuple{tType,tType}
    zspan::Tuple{zType,zType}
    callback::C
    function NLSEProblem(N, D, u0, tspan, zspan; callback=nothing)
        new{typeof(u0), promote_type(map(typeof, tspan)...),
            promote_type(map(typeof, zspan)...),
            typeof(N), typeof(D), typeof(callback)}(
            N, D, u0, tspan, zspan, callback)
    end
end
