module NonlinearOptics
    using Polynomials
    using OrdinaryDiffEq, DiffEqBase
    using Juno
    using RecursiveArrayTools, MuladdMacro, Parameters

    using RecipesBase

    import DiffEqCallbacks: PeriodicCallback
    import DiffEqBase: solve, solve!, init, step!,
                       build_solution, initialize!, isinplace
    import Base: getindex, setindex!

    abstract type AbstractNLSEProblem{uType, tType, zType} <: DEProblem end

    abstract type AbstractNLSESolution{T,N} <: AbstractTimeseriesSolution{T,N} end

    abstract type AbstractNLSEAlgorithm <: DEAlgorithm end

    const c = 2.99792458e8/1e12 # preferred units are THz

    include("interface/fit_util.jl")
    include("interface/types.jl")
    include("interface/solvertypes.jl")
    include("interface/conversions.jl")
    include("interface/interface.jl")
    include("optics/optics_utils.jl")
    include("optics/dispersion.jl")
    include("optics/pulses.jl")
    include("interface/build.jl")
    include("interface/solve.jl")
    include("materials/materials.jl")
    include("plot_recipes.jl")
    include("example_models.jl")

    export NLSEProblem, NLSESolution, SymmetrizedSplitStep
    export prob_bright_soliton, prob_linear_pulse, prob_dispersive_pulse,
            prob_LL
    export solve
    export build_problem

    export Wavelength, Frequency, c
    export OpticalAttr
    export Mode
    export CircularResonator, RacetrackResonator, Waveguide
    export frequency, wavelength, getω, get_beta, get_FSR, get_groupindex,
           get_property, get_label
    export add_mode!
    export ToyModel
    export DynamicNLSE, DynamicLL, DynamicIkeda,
           SteadyStateNLSE, SteadyStateLL, SteadyStateIkeda
    export FT

end # module
