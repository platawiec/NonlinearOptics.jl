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
    corefraction::OpticalAttr
end
Mode(effectiveindex, effectivearea) = Mode(effectiveindex, effectivearea, CoreFraction())

struct ToyMode{betaType, areaType, coreType} <: AbstractMode
    beta::Vector{betaType}
    effectivearea::areaType
    corefraction::coreType
end

abstract type AbstractStructure end

mutable struct Waveguide{T} <: AbstractStructure
    length::T
    orientation::Int
    modes::Vector{Mode}
end
Waveguide(length, orientation) = Waveguide{typeof(length)}(length, orientation, Mode[])

abstract type AbstractResonator <: AbstractStructure end
mutable struct CircularResonator{T} <: AbstractResonator
    radius::T
    modes::Vector{Mode}
end
CircularResonator(radius) = CircularResonator{typeof(radius)}(radius, Mode[])
mutable struct RacetrackResonator{T} <: AbstractResonator
    radius::T
    length::T
    orientation::Int
    modes::Vector{Mode}
end
RacetrackResonator(radius, length, orientation) = RacetrackResonator{typeof(radius)}(radius, length, orientation, Mode[])
