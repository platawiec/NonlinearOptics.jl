abstract type AbstractSource end

struct Frequency{T} <: AbstractSource
    f::T
end
struct Wavelength{T} <: AbstractSource
    λ::T
end

abstract type AbstractMaterial end

abstract type AbstractRamanResponse end
struct RamanCoefficient{T} <: AbstractRamanResponse
    raman::T
end
struct RamanTensor <: AbstractRamanResponse end

abstract type AbstractDielectric end
struct DielectricCoefficient{T} <: AbstractDielectric
    ϵ::T
end
struct DielectricTensor <: AbstractDielectric end

struct Glass{T} <: AbstractMaterial
    nonlinearindex::T
    raman::RamanCoefficient
    ϵ::DielectricCoefficient
end
struct Crystal{T} <: AbstractMaterial
    nonlinearindex::T
    raman::RamanTensor
    ϵ::DielectricTensor
end
Glass(nl_index, r, ϵ) = Glass{typeof(nl_index)}(nl_index, RamanCoefficient(r), DielectricCoefficient(ϵ))
Crystal(nl_index, r, ϵ) = Crystal{typeof(nl_index)}(nl_index, RamanTensor(), DielectricTensor())

abstract type AbstractOpticalAttr end

mutable struct OpticalAttr{T, F} <: AbstractOpticalAttr
    source::Vector{Frequency}
    property::Vector{T}
    label::String
    fit_func::F
    OpticalAttr{T, F}(s, p, label) where {T, F} = fit(new(s, p, label))
    OpticalAttr{T, F}(s, p, label, f) where {T, F} = new(s, p, label, f)
end

function OpticalAttr(source, property, label, fit_func)
    OpticalAttr{eltype(property), typeof(fit_func)}(source, property, label, fit_func)
end
function OpticalAttr(source, property, label)
    T = eltype(property)
    OpticalAttr{T, ScaledFit{T, Poly{T}}}(source, property, label)
end
function OpticalAttr(property::Number, label)
    T = typeof(property)
    OpticalAttr{T, Poly{T}}([Frequency(0.0)], [property], label, Poly([property]))
end

abstract type AbstractMode end
mutable struct Mode <: AbstractMode
    effectiveindex::OpticalAttr
    effectivearea::OpticalAttr
    linearloss::OpticalAttr
    corefraction::OpticalAttr
    coupling::OpticalAttr
end
Mode(n_eff, area_eff) = Mode(n_eff, area_eff,
                             OpticalAttr(0.0, "Linear Loss"),
                             OpticalAttr(1.0, "Core Fraction"),
                             OpticalAttr(1.0, "Coupling"))

struct ToyMode <: AbstractMode
    beta::Vector{OpticalAttr}
    effectivearea::OpticalAttr
    linearloss::OpticalAttr
    corefraction::OpticalAttr
    coupling::OpticalAttr
end
ToyMode(beta, area_eff) = ToyMode(beta, area_eff,
                                  OpticalAttr(0.0, "Linear Loss (1/m)"),
                                  OpticalAttr(1.0, "Core Fraction"),
                                  OpticalAttr(1.0, "Coupling"))

abstract type AbstractStructure end

mutable struct Waveguide{T} <: AbstractStructure
    length::T
    orientation::Int
    material::AbstractMaterial
    modes::Vector{AbstractMode}
end
Waveguide(l, orient, mat) = Waveguide{typeof(l)}(l, orient, mat, AbstractMode[])

abstract type AbstractResonator <: AbstractStructure end
mutable struct CircularResonator{T} <: AbstractResonator
    radius::T
    material::AbstractMaterial
    modes::Vector{AbstractMode}
end
CircularResonator(r, mat) = CircularResonator{typeof(r)}(r, mat, AbstractMode[])
mutable struct RacetrackResonator{T} <: AbstractResonator
    radius::T
    length::T
    orientation::Int
    material::AbstractMaterial
    modes::Vector{AbstractMode}
end
RacetrackResonator(r, l, mat, orient) = RacetrackResonator{typeof(r)}(
                                                    r, l, mat, orient, AbstractMode[])

abstract type AbstractLaser end
mutable struct CWLaser{T1, T2} <: AbstractLaser
    frequency::Frequency
    detuning::T1
    power::T2
end
CWLaser(f, δ, P) = CWLaser{typeof(δ), typeof(P)}(f, δ, P)
mutable struct PulsedLaser <: AbstractLaser
    frequency::Frequency
    pulse_init
end
PulsedLaser(f, avgP, reprate, pulsetime, kwargs...) = PulsedLaser(
                                            f, derive_pulse(avgP, reprate, pulsetime, kwargs...))

abstract type AbstractModel end
mutable struct Model <: AbstractModel
    laser::AbstractLaser
    structure::AbstractStructure
end
mutable struct ToyModel{T} <: AbstractModel
    ω0::T
    FSR::T
    nonlinearcoeff::T
    linearloss::T
    coupling::T
    power_in::T
    length::T
    detuning::T
    betacoeff::Vector{T}
    pulsetime::T
    has_shock::Bool
    has_raman::Bool
end
ToyModel(;ω0=200., FSR=0.1, nonlinearcoeff=1.0, linearloss=0.009, coupling=0.009,
          power_in = 0.755, length=628e-6, detuning=0.0534, pulsetime=0.1,
          betacoeff=[0, 0, -0.05],
          has_shock=false, has_raman=false) = ToyModel(ω0, FSR, nonlinearcoeff, linearloss,
                                              coupling, power_in, length,
                                              detuning, betacoeff, pulsetime,
                                              has_shock, has_raman)
