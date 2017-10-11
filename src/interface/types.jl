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

abstract type AbstractOpticalProperty{T} end

struct GenericOpticalProperty{T} <: AbstractOpticalProperty{T}
    source::Vector{AbstractSource}
    property::Vector{T}
    label::String
    fit_func::ScaledFit{T, Poly{T}}
    GenericOpticalProperty{T}(s, p, label, f) where T = fit(new(s, p, label))
end
function GenericOpticalProperty(source, property, label)
    GenericOpticalProperty{eltype(property)}(source, property, label)
end

mutable struct EffectiveRefractiveIndex{T} <: AbstractOpticalProperty{T}
    source::Vector{AbstractSource}
    effectiveindex::Vector{T}
    fit_func::ScaledFit{T, Poly{T}}
    EffectiveRefractiveIndex{T}(s, p) where T = fit(new(s, p))
end
function EffectiveRefractiveIndex(source, property)
    EffectiveRefractiveIndex{eltype(property)}(source, property)
end

struct CoreFraction{T} <: AbstractOpticalProperty{T}
    source::Vector{AbstractSource}
    corefraction::Vector{T}
    fit_func::ScaledFit{T, Poly{T}}
    CoreFraction{T}(s, p) where T = fit(new(s, p))
end
CoreFraction() = CoreFraction([Wavelength(0.0)], [1.0])
CoreFraction(s, p) = CoreFraction{eltype(p)}(s, p)

struct EffectiveModeArea{T} <: AbstractOpticalProperty{T}
    source::Vector{AbstractSource}
    effectivearea::Vector{T}
    fit_func::ScaledFit{T, Poly{T}}
    EffectiveModeArea{T}(s, p) where T = fit(new(s, p))
end
EffectiveModeArea(s, p) = EffectiveModeaArea{eltype(p)}(s, p)

abstract type AbstractMode end
struct Mode <: AbstractMode
    effectiveindex::EffectiveRefractiveIndex
    effectivearea::EffectiveModeArea
    corefraction::CoreFraction
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
