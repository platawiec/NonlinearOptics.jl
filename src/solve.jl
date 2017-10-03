function solve{algType<:AbstractNLSEAlgorithm}(
    prob::AbstractNLSEProblem,
    alg::algType,
    timeseries=[],
    ts=[], ks=[];
    kwargs...)

    integrator = init(prob, alg, timeseries, ts, prob.zmesh, ks; kwargs...)
    solve!(integrator)
    integrator.sol
end

function init{algType<:AbstractNLSEAlgorithm}(
    prob::AbstractNLSEProblem,
    alg::algType,
    timeseries_init=typeof(prob.u0)[],
    ts_init=eltype(prob.tspan)[], zmesh_init=eltype(prob.zspan)[],
    ks_init=[];
    timeseries_steps = 1,
    saveat = eltype(prob.tspan)[],
    tstops = eltype(prob.tspan)[],
    save_everystep = isempty(saveat),
    save_timeseries = nothing,
    save_start = true,
    dense = save_everystep,
    dt = (prob.tspan[2]-prob.tspan[1])/10000,
    maxiters = 1000000,
    verbose = true,
    progress = false,
    progress_steps = 1000,
    progress_name = "NLSE",
    progress_message = "NLSE",
    kwargs...
    )

    progress ? (prog = Juno.ProgressBar(name=progress_name)) : prog = nothing

    u = deepcopy(prob.u0)
    uType = typeof(u)
    uprev = deepcopy(u)
    t = prob.tspan[1]
    tType = typeof(t)
    tstops = collect(prob.tspan[1]:dt:prob.tspan[2])
    zs = zmesh
    ktilde = 2pi./zs

    ks = ks_init
    timeseries = timeseries_init

    cache = alg_cache(alg, u)

    CacheType = typeof(cache)

    sol = build_solution(prob, alg, tstops, zs, timeseries)

    integrator = NLSEIntegrator{algType, uType, tType, typeof(tstops), typeof(ktilde), CacheType, typeof(prob.N), typeof(prob.D)}(
        N, D, sol, u, uprev, t, dt, tstops, ktilde, alg, cache)

    initialize!(integrator, integrator.cache)

    return integrator
end

function solve!(integrator::NLSEIntegrator)
    @inbounds while !isempty(integrator.tstops)
        integrator.t = pop!(integrator.tstops)
        perform_step!(integrator, integrator.cache)
        push!(integrator.sol.u, integrator.u)
        #push!(integrator., integrator.u)
    end


    integrator.sol.retcode = :Success
    nothing
end
