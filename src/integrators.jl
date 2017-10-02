mutable struct NLSEIntegrator{algType<:AbstractNLSEAlgorithm, uType, tType, CacheType<:DECache}
    sol::NLSESolution
    u::uType
    uprev::uType
    t::tType
    dt::tType
    alg::algType
    cache::CacheType
end

# placeholder for future cache-ing
function initialize!(integrator, cache::SplitStepConstantCache)

end

function perform_step!(integrator, cache::SplitStepConstantCache, repeat_step=false)
    @unpack t, dt, uprev, u = integrator
    @unpack fft_planned!, ifft_planned! = cache

    halfdt = dt/2

    @. u = exp(N(uprev) * halfdt)
    u = fftshift(fft_planned! * u)
    @. u = exp(D(uprev) * dt)
    u = ifft_planned! * ifftshift(u)
    @. u = exp(N(uprev) * halfdt)

    integrator.uprev = copy(u)
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
