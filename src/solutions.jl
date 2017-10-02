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
(sol::NLSESolution)(t, z, deriv::Type=Val{0}; idx=nothing) = sol.interp(t, z, idxs, deriv)
(sol::NLSESolution)(v, t, z, deriv::Type=Val{0}; idx=nothing) = sol.interp(v, t, z, idxs, deriv)

function build_solution(
    prob::AbstractNLSEProblem,
    alg, t, z, u;
    dense=false, k=[], du=[],
    interp = !isempty(du) ? HermiteInterpolation(t, u, du) : LinearInterpolation(t, u),
    retcode = :Default, kwargs...)

    T = eltype(eltype(u))
    N = length(size(prob.u0)..., length(u))

    sol = NLSESolution{T,N,typeof(u),Void,Void,
                       typeof(t),typeof(z), typeof(k),
                       typeof(prob),typeof(alg),typeof(interp)}(u,nothing,nothing,
                       t,k,prob,alg,interp,dense,0,retcode)

    return sol
end
