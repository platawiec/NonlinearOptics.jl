module NonlinearOptics
    using DiffEqBase
    using Juno
    using RecursiveArrayTools, MuladdMacro, Parameters

    using RecipesBase

    import DiffEqBase: solve, solve!, init, step!,
                       build_solution, initialize!, isinplace

    abstract type AbstractNLSEProblem{uType, tType, zType} <: DEProblem end

    abstract type AbstractNLSESolution{T,N} <: AbstractTimeseriesSolution{T,N} end

    abstract type AbstractNLSEAlgorithm <: DEAlgorithm end

    include("problems.jl")
    include("solutions.jl")
    include("algorithms.jl")
    include("integrators.jl")
    include("solve.jl")
    include("premade_problems.jl")
    include("plot_recipes.jl")

    export NLSEProblem

    export NLSESolution

    export AbstractNLSEAlgorithm, SymmetrizedSplitStep

    export prob_bright_soliton, prob_linear_pulse

    export solve

end # module
