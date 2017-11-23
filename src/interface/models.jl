# defines interfaces and traits of Model and ToyModel types
pathlength(w::Waveguide) = w.length
pathlength(res::CircularResonator) = 2pi*res.radius
pathlength(res::RacetrackResonator) = 2pi*res.radius + 2*res.length

pathlength(m::ToyModel) = m.length
pathlength(m::Model) = sum(pathlength, m.structure)

get_ω0(m::ToyModel) = m.ω0
get_ω0(m::Model) = getω(m.laser)

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
    loss = zeros(eltype(ω), length(ω), num_modes(model), num_structures(model))
    for (i, structure) in enumerate(model.structure)
        for (j, mode) in enumerate(structure.modes)
            loss[:, j, i] = mode.linearloss.(ω)
        end
    end
    return loss
end

nonlinearcoeff(model::ToyModel, ω) = model.nonlinearcoeff
function nonlinearcoeff(model::Model, ω)
    #TODO: fix up zero setting
    # pull n2 from structure
    nlcoeff = zeros(eltype(ω), length(ω), num_modes(model), num_structures(model))
    for (i, structure) in enumerate(model.structure)
        nl_n = structure.material.nonlinearindex
        for (j, mode) in enumerate(structure.modes)
            nlcoeff[:, j, i] = ω .* nl_n .* mode.corefraction.(ω) ./ mode.effectivearea.(ω) / c
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
    dispersion = zeros(eltype(ω), length(ω), num_modes(model), num_structures(model))
    for (i, structure) in enumerate(model.structure)
        for (j, mode) in enumerate(structure.modes)
            dispersion[:, j, i] = get_stationarydispersion(mode, model.laser).(ω)
        end
    end
    return dispersion
end

has_raman(m::ToyModel) = m.has_raman
function has_raman(m::Model)
    has_raman = fill(false, num_modes(m))
    for structure in m.structure
        for (i, mode) in enumerate(structure.modes)
            has_raman[i] = mode.has_raman
        end
    end
    return has_raman
end

num_structures(m::Model) = length(m.structure)
num_modes(m::Model) = length(m.structure[1].modes)

function get_structure_idx(m::Model, z)
    structure_pos = cumsum(pathlength.(m.structure))
    idx = count(i->(z>i), structure_pos)+1
    return idx
end
get_structure_idx(m::ToyModel, z) = 1
# placeholder for now
function get_orientation(structure, z)
    return Vec{3}((1, 0, 0)), π/4
end

# need to get material tensor and wait for Tensors.jl to merge with rotation
function raman_tensor(structure, z)
    axis_vec, rot_angle = get_orientation(structure, z)
    #tensor = material.raman_tensor
    tensor = CubicRamanTensor() # placeholder
    #tensor = rotate(tensor, axis_vec, rot_angle)
    return tensor
end
function electronic_tensor(structure, z)
    axis_vec, rot_angle = get_orientation(structure, z)
    #tensor = material.electronic_tensor
    tensor = CubicElectronicTensor(0.1) # TODO placeholders
    #tensor = rotate(tensor, axis_vec, rot_angle)
    return tensor
end

interacts_with_other_mode(mode) = mode.has_interaction
get_paired_mode_idx(structure, mode_idx) = structure.interactions[mode_idx]
get_polarization(mode) = mode.polarization == :TM ? 1 : 2
