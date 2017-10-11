abstract type AbstractSource end

struct Frequency <: AbstractSource
    f::Real
end
struct Wavelength <: AbstractSource
    λ::Real
end

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

abstract type AbstractOpticalProperty end

struct GenericOpticalProperty{T} <: AbstractOpticalProperty
    source::Vector{AbstractSource}
    property::Vector{T}
    fit_func::ScaledFit{T, Poly}
    label::String
end
function GenericOpticalProperty(source, property, label)
    return GenericOpticalProperty{eltype(property)}(
        source, property, get_fit, name
    )
end

struct EffectiveRefractiveIndex{T} <: AbstractOpticalProperty{T}
    source::Vector{AbstractSource}
    effectiveindex::Vector{T}
    fit_func::ScaledFit{T, Poly}
end
EffectiveRefractiveIndex(s, prop) = EffectiveRefractiveIndex{eltype(prop)}(s, prop)

struct CoreFraction <: AbstractOpticalProperty
    source::Vector{AbstractSource}
    corefraction::Vector{T}
    fit_func::ScaledFit{T, Poly}
end
CoreFraction() = CoreFraction([Wavelength(0.0)], [1.0])

struct EffectiveModeArea{T} <: AbstractOpticalProperty{T}
    source::Vector{AbstractSource}
    effectivearea::Vector{T}
    fit_func::ScaledFit{T, Poly}
end

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

struct Waveguide <: AbstractStructure
    length::Real
    orientation::Int
    modes::Vector{Mode}
end
abstract type AbstractResonator <: AbstractStructure end
struct CircularResonator <: AbstractResonator
    radius::Real
    modes::Vector{Mode}
end
struct RacetrackResonator <: AbstractResonator
    radius::Real
    length::Real
    orientation::Int
    modes::Vector{Mode}
end
