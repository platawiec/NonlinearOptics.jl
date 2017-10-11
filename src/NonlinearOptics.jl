module NonlinearOptics
    using Polynomials
    using DiffEqBase
    using Juno
    using RecursiveArrayTools, MuladdMacro, Parameters

    using RecipesBase

    import DiffEqBase: solve, solve!, init, step!,
                       build_solution, initialize!, isinplace

    abstract type AbstractNLSEProblem{uType, tType, zType} <: DEProblem end

    abstract type AbstractNLSESolution{T,N} <: AbstractTimeseriesSolution{T,N} end

    abstract type AbstractNLSEAlgorithm <: DEAlgorithm end

    const c = 2.99792458e8

    include("interface/fit_util.jl")
    include("interface/types.jl")
    include("interface/interface.jl")
    include("diffeqs/problems.jl")
    include("diffeqs/algorithms.jl")
    include("diffeqs/solutions.jl")
    include("diffeqs/algorithms.jl")
    include("diffeqs/integrators.jl")
    include("diffeqs/solve.jl")
    include("diffeqs/premade_problems.jl")
    include("plot_recipes.jl")

    export NLSEProblem

    export NLSESolution

    export AbstractNLSEAlgorithm, SymmetrizedSplitStep

    export prob_bright_soliton, prob_linear_pulse, prob_dispersive_pulse,
            prob_LL

    export solve

end # module
