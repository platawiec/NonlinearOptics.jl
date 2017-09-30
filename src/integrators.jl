struct SymmetrizedSplitStep <: AbstractNLSEAlgorithm end

# placeholder for future cache-ing
function initialize!(integrator)

end

function perform_step!(integrator, repeat_step=false)
    @unpack t, dt, z, dz, uprev, u, f, fft_planned = integrator

    dzhalf = dz/2.0

    @. u = exp(N(uprev) * dzhalf)
    u = fftshift(fft_planned * u)
    @. u = exp(D(uprev) * dz)
    u = ifft_planned * ifftshift(u)
    @. u = exp(N(uprev) * dzhalf)

    integrator.uprev = u

end
