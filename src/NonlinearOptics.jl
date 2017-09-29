module NonlinearOptics
    using DiffEqPDEBase
    using DiffEqBase
    using Juno
    using RecursiveArrayTools, MuladdMacro

    import DiffEqPDEBase: solve

    abstract type AbstractNLSEProblem{uType, tType, zType} <: PDEProblem end

    abstract type AbstractNLSESolution{T,N} <: AbstractTimeseriesSolution{T,N} end

    abstract type AbstractNLSEAlgorithm <: PDEAlgorithm end

    include("problems.jl")
    include("solutions.jl")
    include("algorithms.jl")
    include("integrators.jl")
    include("solve.jl")
    include("premade_problems.jl")

    export NLSEProblem

    export NLSESolution

    export NLSEAlgorithm

    export prob_bright_soliton


end # module
