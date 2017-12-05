struct RamanTensor{T, T1, T2, T3}
    tensor::T
    tau1::T1
    tau2::T2
    raman_fraction::T3
end
struct ElectronicTensor{T, T1}
    tensor::T
    nl_index::T1
end

abstract type AbstractDielectric end
struct DielectricCoefficient <: AbstractDielectric end
struct DielectricTensor <: AbstractDielectric end

struct Material
    electronic::ElectronicTensor
    raman::RamanTensor
end

mutable struct OpticalAttr{TF<:Unitful.Frequency, T, F}
    source::Vector{TF}
    property::Vector{T}
    label::String
    fit_func::F
    OpticalAttr{TF, T, F}(s, p, label) where {TF, T, F} = fit(new(s, p, label))
    OpticalAttr{TF, T, F}(s, p, label, f) where {TF, T, F} = new(s, p, label, f)
end

function OpticalAttr(source, property, label, fit_func)
    OpticalAttr{eltype(source), eltype(property), typeof(fit_func)}(source, property, label, fit_func)
end
function OpticalAttr(source::AbstractVector{TF}, property, label) where {TF<:Unitful.Frequency}
    T = eltype(property)
    OpticalAttr{TF, T, ScaledFit{TF, Poly{T}}}(source, property, label)
end
function OpticalAttr(source::AbstractVector{TL}, property, label) where {TL<:Unitful.Length}
    source_infrequency = frequency.(source)
    TF = eltype(source_infrequency)
    T = eltype(property)
    OpticalAttr{TF, T, ScaledFit{TF, Poly{T}}}(source_infrequency, property, label)
end
function OpticalAttr(property::Number, label)
    TF = typeof(1.0THz2π)
    T = typeof(property)
    OpticalAttr{TF, T, Poly{T}}([TF(0.0)], [property], label, Poly([property]))
end

abstract type AbstractMode end
mutable struct Mode <: AbstractMode
    effectiveindex::OpticalAttr
    effectivearea::OpticalAttr
    linearloss::OpticalAttr
    corefraction::OpticalAttr
    coupling::OpticalAttr
    polarization::Symbol
    has_shock::Bool
    has_raman::Bool
    has_interaction::Bool
end
Mode(n_eff, area_eff) = Mode(n_eff, area_eff,
                             OpticalAttr(0.0, "Linear Loss"),
                             OpticalAttr(1.0, "Core Fraction"),
                             OpticalAttr(1.0, "Coupling"),
                             :TM,
                             false,
                             false,
                             false,
                             )
Mode(n_eff, area_eff, lin_loss, corefrac, coupling, pol, shock, raman) = Mode(
     n_eff, area_eff, lin_loss, corefrac, coupling, pol, shock, raman, false)

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

abstract type AbstractLaser end
mutable struct CWLaser{TF<:Unitful.Frequency, T1, TP<:Unitful.Power} <: AbstractLaser
    frequency::typeof(1.0THz2π)
    detuning::T1
    power::TP
end
CWLaser(f, δ, P) = CWLaser{typeof(δ), typeof(P)}(f, δ, P)
mutable struct PulsedLaser{TF<:Unitful.Frequency,
                           TP<:Unitful.Power,
                           TT<:Unitful.Time} <: AbstractLaser
    frequency::TF
    power_in::TP
    pulsetime::TT
end
