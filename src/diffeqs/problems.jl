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
    zmesh::Vector{zType}
    callback::C
end

function NLSEProblem(N, D, u0, tspan, zmesh; callback=nothing, kwargs...)
    u0_vec = u0.(zmesh)
    NLSEProblem{typeof(u0_vec), promote_type(map(typeof, tspan)...),
                eltype(zmesh),
                typeof(N), typeof(D), typeof(callback)}(
                N, D, u0_vec, tspan, zmesh, callback)
end

isinplace(prob::AbstractNLSEProblem{uType, tType, zType}) where {uType, tType, zType} = false
