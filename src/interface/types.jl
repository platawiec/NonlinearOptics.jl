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

struct Material <: AbstractMaterial
    raman::AbstractRamanResponse
    ϵ::AbstractDielectric
end

abstract type AbstractOpticalProperty end

struct EffectiveRefractiveIndex <: AbstractOpticalProperty
    light::Vector{Wavelength}
    effectiveindex::Vector{Real}
end
struct CoreFraction <: AbstractOpticalProperty
    light::Vector{AbstractLight}
    corefraction::Vector{Real}
end
CoreFraction() = CoreFraction([Wavelength(0.0)], [1.0])

struct EffectiveModeArea <: AbstractOpticalProperty
    light::Vector{AbstractLight}
    effectivearea::Vector{Real}
end

struct Mode
    effectiveindex::EffectiveRefractiveIndex
    effectivearea::EffectiveModeArea
    corefraction::CoreFraction
end
Mode(effectiveindex, effectivearea) = Mode(effectiveindex, effectivearea, CoreFraction())

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
