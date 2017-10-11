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
    include("interface/conversions.jl")
    include("interface/interface.jl")
    include("diffeqs/problems.jl")
    include("diffeqs/algorithms.jl")
    include("diffeqs/solutions.jl")
    include("diffeqs/algorithms.jl")
    include("diffeqs/integrators.jl")
    include("diffeqs/solve.jl")
    include("diffeqs/premade_problems.jl")
    include("plot_recipes.jl")
    include("example_models.jl")

    export NLSEProblem, NLSESolution, SymmetrizedSplitStep
    export prob_bright_soliton, prob_linear_pulse, prob_dispersive_pulse,
            prob_LL
    export solve

    export Wavelength, Frequency, c
    export GenericOpticalProperty, EffectiveModeArea, EffectiveRefractiveIndex,
           CoreFraction
    export Mode
    export CircularResonator, RacetrackResonator, Waveguide
    export frequency, wavelength, getω, get_beta, get_FSR, get_groupindex,
           get_property, get_label

end # module
