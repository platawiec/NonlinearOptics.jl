abstract type AbstractSource end

struct Frequency{T} <: AbstractSource
    f::T
end
Frequency(f) = Frequency{typeof(f)}(f)
struct Wavelength{T} <: AbstractSource
    λ::T
end
Wavelength(λ) = Wavelength{typeof(λ)}(λ)

abstract type AbstractMaterial end

abstract type AbstractRamanResponse end

struct RamanCoefficient <: AbstractRamanResponse end
struct RamanTensor <: AbstractRamanResponse end

abstract type AbstractDielectric end

struct DielectricCoefficient <: AbstractDielectric end
struct DielectricTensor <: AbstractDielectric end

struct Material <: AbstractMaterial
    raman::AbstractRamanResponse
    ϵ::AbstractDielectric
end

abstract type AbstractOpticalAttr{T} end

mutable struct OpticalAttr{T} <: AbstractOpticalAttr{T}
    source::Vector{AbstractSource}
    property::Vector{T}
    label::String
    fit_func::ScaledFit{T, Poly{T}}
    OpticalAttr{T}(s, p, label) where T = fit(new(s, p, label))
end
function OpticalAttr(source, property, label)
    OpticalAttr{eltype(property)}(source, property, label)
end

abstract type AbstractMode end
mutable struct Mode <: AbstractMode
    effectiveindex::OpticalAttr
    effectivearea::OpticalAttr
    corefraction::OpticalAttr
end
Mode(effectiveindex, effectivearea) = Mode(effectiveindex, effectivearea, CoreFraction())

struct ToyMode{betaType, areaType, coreType} <: AbstractMode
    beta::Vector{betaType}
    effectivearea::areaType
    corefraction::coreType
end

abstract type AbstractStructure end

mutable struct Waveguide <: AbstractStructure
    length::Real
    orientation::Int
    modes::Vector{Mode}
end
abstract type AbstractResonator <: AbstractStructure end
mutable struct CircularResonator <: AbstractResonator
    radius::Real
    modes::Vector{Mode}
end
mutable struct RacetrackResonator <: AbstractResonator
    radius::Real
    length::Real
    orientation::Int
    modes::Vector{Mode}
end
