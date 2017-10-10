abstract type AbstractLight end

struct Frequency <: AbstractLight
    f::Real
end
struct Wavelength <: AbstractLight
    λ::Real
end

abstract type AbstractMaterial end

abstract type AbstractRamanResponse end

struct RamanCoefficient <: AbstractRamanResponse end
struct RamanTensor <: AbstractRamanResponse end

abstract type AbstractDielectric end

struct DielectricCoefficient <: AbstractDielectric end
struct DielectricTensor <: AbstractDielectric end

abstract type AbstractOpticalProperty end

struct EffectiveRefractiveIndex <: AbstractOpticalProperty
    light::Vector{AbstractLight}
    effectiveindex::Vector{Real}
end
struct CoreFraction <: AbstractOpticalProperty
    light::Vector{AbstractLight}
    corefraction::Vector{Real}
end
struct EffectiveModeArea <: AbstractOpticalProperty
    light::Vector{AbstractLight}
    effectivearea::Vector{Real}
end

abstract type AbstractStructure end

struct Waveguide <: AbstractStructure
    length::Real
    orientation::Int
    effectiveindex::EffectiveRefractiveIndex
    effectivearea::EffectiveModeArea
    corefraction::CoreFraction
end
abstract type AbstractResonator <: AbstractStructure end
struct CircularResonator <: AbstractResonator
    radius::Real
    effectiveindex::EffectiveRefractiveIndex
    effectivearea::EffectiveModeArea
    corefraction::CoreFraction
end
struct RacetrackResonator <: AbstractResonator
    radius::Real
    length::Real
    orientation::Int
    effectiveindex::EffectiveRefractiveIndex
    effectivearea::EffectiveModeArea
    corefraction::CoreFraction
end
