function derive_pulse(average_power, rep_rate, pulse_time, pulse_type::Symbol=:gaussian, chirp=0)
    pulse_energy = sqrt(average_power / rep_rate / pulse_time)
    #u0(x) = (1+0im) * pulse_energy * exp(x^2/pulse_time^2)
    u0(x) = (1+0im) * pulse_energy * sech(x/pulse_time)
    return u0
end

function derive_pulse(peak_power, pulse_FWHM)
    sqrtpower = sqrt(peak_power)
    pulse_0 = pulse_FWHM/(2*log(1+sqrt(2)))
    u0(x) = (1+0im)*sqrtpower*sech(x/pulse_0)
end
