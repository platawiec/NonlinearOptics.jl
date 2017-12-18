# Dispatch types to call different problem builders
abstract type AbstractExperiment end
abstract type AbstractDynamicExperiment <: AbstractExperiment end
abstract type AbstractStochasticExperiment <: AbstractExperiment end
abstract type AbstractSteadyStateExperiment <: AbstractExperiment end
struct DynamicLL <: AbstractDynamicExperiment end
struct DynamicNLSE <: AbstractDynamicExperiment end
struct DynamicIkeda <: AbstractDynamicExperiment end

struct StochasticNLSE <: AbstractStochasticExperiment end

struct SteadyStateLL <: AbstractSteadyStateExperiment end
struct SteadyStateNLSE <: AbstractSteadyStateExperiment end
struct SteadyStateIkeda <: AbstractSteadyStateExperiment end

# Abstraction layers over DiffEq solution/problem types
abstract type AbstractNLOProblem end
mutable struct DynamicNLSEProblem{fType, f0Type, fftType, ifftType, meshType, opType} <: AbstractNLOProblem
    prob::ODEProblem
    model::AbstractModel
    ω::fType
    ω0::f0Type
    tmesh::meshType
    planned_fft::fftType
    planned_ifft::ifftType
    Doperator::opType
end
mutable struct StochasticNLSEProblem{fType, f0Type, fftType, ifftType, meshType, opType} <: AbstractNLOProblem
    prob::SDEProblem
    model::AbstractModel
    ω::fType
    ω0::f0Type
    tmesh::meshType
    planned_fft::fftType
    planned_ifft::ifftType
    Doperator::opType
end
mutable struct DynamicLLProblem{fType, f0Type, fftType, ifftType, meshType, opType} <: AbstractNLOProblem
    prob::ODEProblem
    model::AbstractModel
    ω::fType
    ω0::f0Type
    tmesh::meshType
    planned_fft::fftType
    planned_ifft::ifftType
    Doperator::opType
end
mutable struct DynamicIkedaProblem{fType, f0Type, fftType, ifftType, meshType, opType, cbType} <: AbstractNLOProblem
    prob::ODEProblem
    model::AbstractModel
    ω::fType
    ω0::f0Type
    tmesh::meshType
    planned_fft::fftType
    planned_ifft::ifftType
    Doperator::opType
    ikeda_callback::cbType
end

abstract type AbstractNLOSolution end
mutable struct DynamicNLOSolution <: AbstractNLOSolution
    sol::DESolution
    prob::AbstractNLOProblem
end

# calling solution returns time-domain solution
function (sol::DynamicNLOSolution)(z, mode_idx::Int=1)
    structure_idx, ls = currentstructure(sol.prob.model, z)
    return sol.prob.planned_fft * (sol.sol(z/m)[:, mode_idx, structure_idx].*exp.(view(sol.prob.Doperator, :, mode_idx, structure_idx) * z))
end
function FT(sol::DynamicNLOSolution, z, mode_idx::Int=1)
    dt = sol.prob.tmesh[2]-sol.prob.tmesh[1]
    structure_idx, ls = currentstructure(sol.prob.model, z)
    return fftshift(view(sol.sol(z), :, mode_idx, structure_idx) .* exp.(view(sol.prob.Doperator, :, mode_idx, structure_idx) * z)) / dt
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
