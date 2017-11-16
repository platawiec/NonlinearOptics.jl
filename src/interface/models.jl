# defines interfaces and traits of Model and ToyModel types
pathlength(w::Waveguide) = w.length
pathlength(res::CircularResonator) = 2pi*res.radius
pathlength(res::RacetrackResonator) = 2pi*res.radius + 2*res.length

pathlength(m::ToyModel) = m.length
pathlength(m::Model) = sum(pathlength, m.structure)

get_ω0(m::ToyModel) = m.ω0
get_ω0(m::Model) = getω(m.source)

derive_pulse(m::Model, t) = derivepulse(m.source, t)
function derive_pulse(m::ToyModel, t)
    sqrtpower = sqrt(m.power_in)
    pulse_0 = m.pulsetime/(2*log(1+sqrt(2)))
    return (1+0im)*sqrtpower*sech(t/pulse_0)
end
