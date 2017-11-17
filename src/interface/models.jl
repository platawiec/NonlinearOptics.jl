# defines interfaces and traits of Model and ToyModel types
pathlength(w::Waveguide) = w.length
pathlength(res::CircularResonator) = 2pi*res.radius
pathlength(res::RacetrackResonator) = 2pi*res.radius + 2*res.length

pathlength(m::ToyModel) = m.length
pathlength(m::Model) = sum(pathlength, m.structure)

get_ω0(m::ToyModel) = m.ω0
get_ω0(m::Model) = getω(m.laser)

derive_pulse(m::Model, t) = derive_pulse(m.laser, t)
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
    loss = zeros(ω)
    for structure in model.structure
        for mode in structure.modes
            loss = mode.linearloss.(ω)
        end
    end
    return loss
end

nonlinearcoeff(model::ToyModel, ω) = model.nonlinearcoeff
function nonlinearcoeff(model::Model, ω)
    #TODO: fix up zero setting
    # pull n2 from structure
    nlcoeff = zero(ω)
    for structure in model.structure
        nl_n = structure.material.nonlinearindex
        for mode in structure.modes
            nlcoeff = nl_n * mode.corefraction.(ω) ./ mode.effectivearea.(ω)
        end
    end
    return nlcoeff
end

function stationarydispersion(model::ToyModel, ω)
    beta = model.betacoeff
    betas = beta .* [1/factorial(n) for n=0:(length(beta)-1)]
    dispersion = Poly(betas).(ω)
    return dispersion
end
function stationarydispersion(model::Model, ω)
    #TODO: finish
    dispersion = zeros(ω)
    for structure in model.structure
        for mode in structure.modes
            dispersion = get_stationarydispersion(mode, model.laser).(ω)
        end
    end
    return dispersion
end

has_raman(m::ToyModel) = m.has_raman
function has_raman(m::Model)
    has_raman = false
    for structure in m.structure
        for mode in structure.modes
            has_raman = mode.has_raman
        end
    end
    return has_raman
end
