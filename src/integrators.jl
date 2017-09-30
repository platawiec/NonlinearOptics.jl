mutable struct NLSEIntegrator{algType<:AbstractNLSEAlgorithm}
    sol::SolType
    alg::algType
    cache::CacheType
end

# placeholder for future cache-ing
function initialize!(integrator)

end

function perform_step!(integrator, repeat_step=false)
    @unpack t, dt, z, dz, dzhalf, uprev, u, f, fft_planned, ifft_planned = integrator

    @. u = exp(N(uprev) * dzhalf)
    u = fftshift(fft_planned * u)
    @. u = exp(D(uprev) * dz)
    u = ifft_planned * ifftshift(u)
    @. u = exp(N(uprev) * dzhalf)

    integrator.uprev = u

end
