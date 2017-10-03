mutable struct NLSEIntegrator{algType<:AbstractNLSEAlgorithm, uType, tType, tstopsType, ktildeType, CacheType<:DECache, F1, F2}
    N::F1
    D::F2
    sol::NLSESolution
    u::uType
    uprev::uType
    t::tType
    dt::tType
    tstops::tstopsType
    ktilde::ktildeType
    alg::algType
    cache::CacheType
end

# placeholder for future cache-ing
function initialize!(integrator, cache::SplitStepConstantCache)

end

function perform_step!(integrator, cache::SplitStepConstantCache, repeat_step=false)
    @unpack t, dt, ktilde, uprev, u, N, D = integrator
    @unpack planned_fft!, planned_ifft! = cache

    halfdt = dt/2
    @. u = exp(N(t, uprev) * halfdt)
    u = fftshift(planned_fft! * u)
    @. u = exp(D(ktilde, u) * dt)
    u = planned_ifft! * ifftshift(u)
    @. u = exp(N(t, uprev) * halfdt)

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
