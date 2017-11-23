module NonlinearOptics
    using Polynomials
    using OrdinaryDiffEq, DiffEqBase, Tensors
    using Juno
    using RecursiveArrayTools, MuladdMacro, Parameters

    using RecipesBase

    import DiffEqCallbacks: PeriodicCallback
    import Base: getindex, setindex!
    import DiffEqBase: solve

    abstract type AbstractNLSEProblem{uType, tType, zType} <: DEProblem end

    abstract type AbstractNLSESolution{T,N} <: AbstractTimeseriesSolution{T,N} end

    abstract type AbstractNLSEAlgorithm <: DEAlgorithm end

    const c = 2.99792458e8/1e12 # preferred units are THz
    const ħ = 1.0545718e-34*1e12 # J⋅ps
    const ħ_over_kBT = 0.0254607837 # T = 300 K, ps

    include("interface/fit_util.jl")
    include("interface/types.jl")
    include("interface/models.jl")
    include("interface/solvertypes.jl")
    include("interface/conversions.jl")
    include("interface/interface.jl")
    include("optics/utils.jl")
    include("optics/dispersion.jl")
    include("optics/pulses.jl")
    include("interface/build.jl")
    include("interface/solve.jl")
    include("materials/materials.jl")
    include("plot_recipes.jl")

    export NLSEProblem, NLSESolution, SymmetrizedSplitStep
    export prob_bright_soliton, prob_linear_pulse, prob_dispersive_pulse,
            prob_LL
    export solve
    export build_problem

    export Wavelength, Frequency, c
    export OpticalAttr
    export Mode
    export CircularResonator, RacetrackResonator, Waveguide
    export PulsedLaser, CWLaser
    export frequency, wavelength, getω, get_beta, get_FSR, get_groupindex,
           get_property, get_label
    export add_mode!, add_interaction!
    export Model, ToyModel
    export DynamicNLSE, DynamicLL, DynamicIkeda,
           SteadyStateNLSE, SteadyStateLL, SteadyStateIkeda
    export FT
    export Glass, Silicon, Diamond

end # module
