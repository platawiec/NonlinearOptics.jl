module NonlinearOptics
    using Polynomials
    using OrdinaryDiffEq, DiffEqBase, StochasticDiffEq, Tensors
    using RecursiveArrayTools, Parameters

    import Unitful
    import Unitful: ps, m, J, THz, °, THz2π, W
    import Unitful: c0, ħ, k

    using RecipesBase

    import DiffEqCallbacks: PeriodicCallback
    import Base: getindex, setindex!
    import DiffEqBase: solve

    abstract type AbstractNLSEProblem{uType, tType, zType} <: DEProblem end

    abstract type AbstractNLSESolution{T,N} <: AbstractTimeseriesSolution{T,N} end

    abstract type AbstractNLSEAlgorithm <: DEAlgorithm end

    const ħ_over_kBT = ħ/(k*300(Unitful.K)) # T = 300 K, ps

    Unitful.preferunits(ps)

    include("interface/fit_util.jl")
    include("interface/types.jl")
    include("waveguides/waveguides.jl")
    include("interface/models.jl")
    include("interface/solvertypes.jl")
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

    export wavelength, frequency, c0
    export OpticalAttr
    export Mode, Waveguide
    export PulsedLaser, CWLaser
    export frequency, wavelength, getω, get_beta, get_FSR, get_groupindex,
           get_property, get_label
    export add_mode!, add_interaction!
    export add_straight!, add_turn!
    export currentmedium, lastmedium, newmedium
    export Model, ToyModel
    export DynamicNLSE, DynamicLL, DynamicIkeda,
           SteadyStateNLSE, SteadyStateLL, SteadyStateIkeda,
           StochasticNLSE
    export FT
    export Glass, Silicon, Diamond

end # module
