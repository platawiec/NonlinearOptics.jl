struct NLSEIntegrator end

# placeholder for future cache-ing
function initialize!(integrator)

end

function perform_step!(integrator, repeat_step=false)
    @unpack t, dt, z, dz, uprev, u, f, fft_planned = integrator

    pulse2 = abs2.(uprev)
    dpulse2 = [0 diff(pulse2)/dt]
    @. nl_mod = exp(1im * gamma * (pulse2 + 1im * tau_shock * dpulse2)*dz)

    pulse = fftshift(plan * pulse)
    @. pulse = dispersion * pulse
    pulse = iplan * ifftshift(pulse)
    @. pulse = nl_mod * pulse
    pulse = fftshift(plan * pulse)
    @. pulse = dispersion * pulse
    pulse = iplan * ifftshift(pulse)
    @. u = abs2(uprev) + dt*integrator.fsalfirst
    f(t+dt, u, integrator.fsallast)

    integrator.u = u
end
