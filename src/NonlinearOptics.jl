module NonlinearOptics
    using DiffEqBase
    using Juno
    using RecursiveArrayTools, MuladdMacro, Parameters

    import DiffEqBase: solve, isinplace

    abstract type AbstractNLSEProblem{uType, tType, zType} <: DEProblem end

    abstract type AbstractNLSESolution{T,N} <: AbstractTimeseriesSolution{T,N} end

    abstract type AbstractNLSEAlgorithm <: DEAlgorithm end

    include("problems.jl")
    include("solutions.jl")
    include("algorithms.jl")
    include("integrators.jl")
    include("solve.jl")
    include("premade_problems.jl")

    export NLSEProblem

    export NLSESolution

    export SymmetrizedSplitStep

    export prob_bright_soliton


end # module
