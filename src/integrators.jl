mutable struct NLSEIntegrator{algType<:AbstractNLSEAlgorithm, CacheType<:DECache}
    sol::NLSESolution
    u::uType
    alg::algType
    cache::CacheType
end

# placeholder for future cache-ing
function initialize!(integrator)

end

function perform_step!(integrator, repeat_step=false)
    @unpack t, dt, z, dz, dzhalf, uprev, u, f, fft_planned!, ifft_planned! = integrator

    @. u = exp(N(uprev) * dzhalf)
    u = fftshift(fft_planned! * u)
    @. u = exp(D(uprev) * dz)
    u = ifft_planned! * ifftshift(u)
    @. u = exp(N(uprev) * dzhalf)

    integrator.u = u

end

function loopheader!(integrator)
end

function loopfooter!(integrator)
    integrator.tprev = integrator.t
    integrator.t = integrator.t + integrator.dt
    if !(typeof(integrator.prog)<:Void) && integrator.opts.progress && integrator.opts.progress_steps==0
        Juno.msg(integrator.prog, integrator.opts.progress_message(integrator.dt, integrator.t, integrator.u))
        Juno.progress(integrator.prog, integrator.t/integrator.sol.prob.tspan[2])
    end
end

function handle_tstop!(integrator)
end

function postamble!(integrator)
end
