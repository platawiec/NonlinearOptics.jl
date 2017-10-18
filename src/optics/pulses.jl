function derive_pulse(average_power, rep_rate, pulse_time, pulse_type::Symbol=:gaussian, chirp=0)
    pulse_energy = sqrt(average_power / rep_rate / pulse_time)
    if pulse_type == :gaussian
        u0(x) = (1+0im) * pulse_energy * exp(x^2/pulse_time^2)
    elseif pulse_type == :secant
        u0(x) = (1+0im) * pulse_energy * sech(x/pulse_time)
    end
    u0
end
