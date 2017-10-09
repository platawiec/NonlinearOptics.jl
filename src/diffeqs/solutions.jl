mutable struct NLSESolution{T,N,uType,uType2,DType,tType,zType,rateType,P,A,IType} <: AbstractNLSESolution{T,N}
  u::uType
  u_analytic::uType2
  errors::DType
  t::tType
  z::zType
  k::rateType
  prob::P
  alg::A
  interp::IType
  dense::Bool
  tslocation::Int
  retcode::Symbol
end
(sol::NLSESolution)(t, deriv::Type=Val{0}; idxs=nothing) = sol.interp(t, idxs, deriv)
(sol::NLSESolution)(v, t, deriv::Type=Val{0}; idxs=nothing) = sol.interp(v, t, idxs, deriv)

function build_solution(
    prob::AbstractNLSEProblem,
    alg, t, z, u;
    dense=false, k=[], du=[],
    interp = !isempty(du) ? DiffEqBase.HermiteInterpolation(t, u, du) : DiffEqBase.LinearInterpolation(t, u),
    retcode = :Default, kwargs...)

    T = eltype(eltype(u))
    N = length(prob.u0)

    sol = NLSESolution{T,N,typeof(u),Void,Void,
                       typeof(t),typeof(z), typeof(k),
                       typeof(prob),typeof(alg),typeof(interp)}(u,nothing,nothing,
                       t,z,k,prob,alg,interp,dense,0,retcode)

    return sol
end
