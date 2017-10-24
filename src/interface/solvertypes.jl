# Dispatch types to call different problem builders
abstract type AbstractExperiment end
abstract type AbstractDynamicExperiment end
abstract type AbstractSteadyStateExperiment end
struct DynamicLL <: AbstractDynamicExperiment end
struct DynamicNLSE <: AbstractDynamicExperiment end
struct DynamicIkeda <: AbstractDynamicExperiment end

struct SteadyStateLL <: AbstractSteadyStateExperiment end
struct SteadyStateNLSE <: AbstractSteadyStateExperiment end
struct SteadyStateIkeda <: AbstractSteadyStateExperiment end

# Abstraction layers over DiffEq solution/problem types
abstract type AbstractNLOProblem end
mutable struct DynamicNLSEProblem{fType, f0Type, fftType, ifftType, meshType, opType} <: AbstractNLOProblem
    prob::ODEProblem
    ω::fType
    ω0::f0Type
    tmesh::meshType
    planned_fft::fftType
    planned_ifft::ifftType
    Doperator::opType
end
mutable struct DynamicLLProblem{fType, f0Type, fftType, ifftType, meshType, opType} <: AbstractNLOProblem
    prob::ODEProblem
    ω::fType
    ω0::f0Type
    tmesh::meshType
    planned_fft::fftType
    planned_ifft::ifftType
    Doperator::opType
end

abstract type AbstractNLOSolution end
mutable struct DynamicNLSESolution <: AbstractNLOSolution
    sol::ODESolution
    prob::DynamicNLSEProblem
end
mutable struct DynamicLLSolution <: AbstractNLOSolution
    sol::ODESolution
    prob::DynamicLLProblem
end


# calling solution returns time-domain solution
function (sol::DynamicNLSESolution)(z)
    sol.prob.planned_fft * (sol.sol(z).*exp(sol.prob.Doperator * z))
end
function FT(sol::DynamicNLSESolution, z)
    dt = sol.prob.tmesh[2]-sol.prob.tmesh[1]
    fftshift(sol.sol(z).*exp(sol.prob.Doperator * z)) / dt
end
function (sol::DynamicLLSolution)(z)
    sol.prob.planned_fft * (sol.sol(z).*exp(sol.prob.Doperator * z))
end
function FT(sol::DynamicLLSolution, z)
    dt = sol.prob.tmesh[2]-sol.prob.tmesh[1]
    fftshift(sol.sol(z).*exp(sol.prob.Doperator * z)) / dt
end

##Abstract interfacei forwards to DiffEq sol field
Base.getindex(A::AbstractNLOSolution,i::Int) = getindex(A.sol, i)
Base.getindex(A::AbstractNLOSolution,I::Vararg{Int, N}) where {N} = getindeX(A.sol, I)
Base.setindex!(A::AbstractNLOSolution, v, i::Int) = setindex!(A.sol, v, i)
Base.setindex!(A::AbstractNLOSolution, v, I::Vararg{Int, N}) where {N} = setindex!(A.sol, v, I)
size(A::AbstractNLOSolution) = size(A.sol)
Base.summary(A::AbstractNLOSolution) = summary(A.sol)
Base.show(io::IO, A::AbstractNLOSolution) = show(io, A.sol)
Base.show(io::IO, m::MIME"text/plain", A::AbstractNLOSolution) = show(io, m, A.sol)
tuples(sol::AbstractNLOSolution) = tuples(sol.sol)
start(sol::AbstractNLOSolution) = start(sol.sol)
next(sol::AbstractNLOSolution,state) = next(sol.sol, state)
done(sol::AbstractNLOSolution,state) = done(sol.sol, state)
