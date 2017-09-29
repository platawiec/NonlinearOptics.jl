"""
𝚤∂ψ/∂t = -1/2∂²ψ/∂x² + κ|ψ|²ψ
"""
mutable struct NLSEProblem{uType, tType, zType} <: AbstractNLSEProblem{uType, tType, zType}
    u0::uType
    tspan::Tuple{tType,tType}
    zspan::Tuple{zType,zType}
end
