abstract type AbstractModel end
mutable struct Model <: AbstractModel
    laser::AbstractLaser
    waveguide::Waveguide
end
mutable struct ToyModel{TF, TG, TL, TPow, TLen, TT} <: AbstractModel
    ω0::TF
    FSR::TF
    nonlinearcoeff::TG
    linearloss::TL
    coupling::Float64
    power_in::TPow
    length::TLen
    detuning::Float64
    betacoeff::Vector
    pulsetime::TT
    has_shock::Bool
    has_raman::Bool
end
ToyModel(;ω0=200.0THz, FSR=0.1THz, nonlinearcoeff=1.0/W/m, linearloss=0.009/m,
          coupling=0.009, power_in = 0.755W, length=(628e-6)m, detuning=0.0534,
          pulsetime=0.1ps, betacoeff=[0/m, 0ps/m, -0.05ps^2/m],
          has_shock=false, has_raman=false) = ToyModel(ω0, FSR, nonlinearcoeff, linearloss,
                                              coupling, power_in, length,
                                              detuning, betacoeff, pulsetime,
                                              has_shock, has_raman)

@inline currentstructure(model::Model, z) = currentstructure(model.waveguide, z)
@inline pathlength(m::Model) = pathlength(m.waveguide)
@inline pathlength(m::ToyModel) = m.length

function derive_pulse(m::Model, t)
    #TODO: Generic injection of pulse to arbitrary mode
    pulse = derive_pulse(m.laser, t)
    full_pulse = zeros(eltype(pulse), length(t), num_modes(m))
    full_pulse[:, 1] = pulse
    return full_pulse
end
function derive_pulse(m::ToyModel, t)
    sqrtpower = sqrt(m.power_in)
    pulse_0 = m.pulsetime/(2*log(1+sqrt(2)))
    return (1+0im)*sqrtpower*sech.(t/pulse_0)
end
function derive_pulse(laser::PulsedLaser, t)
    sqrtpower = sqrt(laser.power_in)
    pulse_0 = laser.pulsetime/(2*log(1+sqrt(2)))
    return (1+0im)*sqrtpower*sech.(t/pulse_0)
end

linearloss(model::ToyModel, ω) = model.linearloss
function linearloss(model::Model, ω)
    firstmode = first(first(model.waveguide.structures).medium.modes)
    TLoss = eltype(firstmode.linearloss.(ω))
    loss = zeros(TLoss, length(ω), num_modes(model), num_structures(model))
    for (i, structure) in enumerate(model.waveguide.structures)
        for (j, mode) in enumerate(structure.medium.modes)
            loss[:, j, i] = mode.linearloss.(ω)
        end
    end
    return loss
end

nonlinearcoeff(model::ToyModel, ω) = model.nonlinearcoeff
function nonlinearcoeff(model::Model, ω)
    firstmode = first(first(model.waveguide.structures).medium.modes)
    nl_n = first(model.waveguide.structures).medium.material.electronic.nl_index
    TNL = eltype(ω .* nl_n .* firstmode.corefraction(first(ω)) / firstmode.effectivearea(first(ω)) / c0)
    # pull n2 from structure
    nlcoeff = zeros(TNL, length(ω), num_modes(model), num_structures(model))
    for (i, structure) in enumerate(model.waveguide.structures)
        nl_n = structure.medium.material.electronic.nl_index
        for (j, mode) in enumerate(structure.medium.modes)
            nlcoeff[:, j, i] = ω .* nl_n .* mode.corefraction.(ω) ./ mode.effectivearea.(ω) / c0
        end
    end
    return nlcoeff
end

function stationarydispersion(model::ToyModel, ω)
    beta = model.betacoeff
    betas = beta .* [1/factorial(n) for n=0:(length(beta)-1)]
    dispersion = Poly(betas).(ω - model.ω0)
    return dispersion
end
function stationarydispersion(model::Model, ω)
    dispersion = zeros(eltype(ω / c0), length(ω), num_modes(model), num_structures(model))
    for (i, structure) in enumerate(model.waveguide.structures)
        for (j, mode) in enumerate(structure.medium.modes)
            dispersion[:, j, i] = get_stationarydispersion(mode, model.laser).(ω)
        end
    end
    return dispersion
end

has_raman(m::ToyModel) = m.has_raman
function has_raman(model::Model)
    has_raman = fill(false, num_modes(model))
    for structure in model.waveguide.structures
        for (i, mode) in enumerate(structure.medium.modes)
            has_raman[i] = mode.has_raman
        end
    end
    return has_raman
end

@inline num_structures(m::Model) = num_structures(m.waveguide)
@inline num_modes(m::Model) = num_modes(m.waveguide)

# need to wait for Tensors.jl to merge with rotation
function raman_tensor(structure, z)
    axis_vec, rot_angle = orientation(structure, z)
    #tensor = material.raman_tensor
    tensor = structure.medium.material.raman.tensor
    #tensor = rotate(tensor, axis_vec, rot_angle)
    return tensor
end
function electronic_tensor(structure, z)
    axis_vec, rot_angle = orientation(structure, z)
    #tensor = material.electronic_tensor
    tensor = structure.medium.material.electronic.tensor
    #tensor = rotate(tensor, axis_vec, rot_angle)
    return tensor
end

@inline function raman_fraction(structure)
    return structure.medium.material.raman.raman_fraction
end

@inline interacts_with_other_mode(mode) = mode.has_interaction
@inline get_paired_mode_idx(med::Medium, mode_idx) = med.interactions[mode_idx]
@inline get_polarization(mode) = mode.polarization == :TM ? 1 : 2
