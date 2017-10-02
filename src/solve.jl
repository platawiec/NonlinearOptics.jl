function solve{algType<:AbstractNLSEAlgorithm}(
    prob::AbstractNLSEProblem,
    alg::algType,
    timeseries=[],
    ts=[], zs=[], ks=[];
    kwargs...)

    integrator = init(prob, alg, timeseries, ts, zs, ks; kwargs...)
    solve!(integrator)
    integrator.sol
end

function init{algType<:AbstractNLSEAlgorithm}(
    prob::AbstractNLSEProblem,
    alg::algType,
    timeseries_init=typeof(prob.u0)[],
    ts_init=eltype(prob.tspan)[], zs_init=eltype(prob.zspan)[],
    ks_init=[];
    timeseries_steps = 1,
    saveat = eltype(prob.tspan)[],
    tstops = eltype(prob.tspan)[],
    save_everystep = isempty(saveat),
    save_timeseries = nothing,
    save_start = true,
    dense = save_everystep,
    dt = (p.tspan[2]-p.tspan[1])/10000,
    maxiters = 1000000,
    verbose = true,
    progress = false,
    progress_steps = 1000,
    progress_name = "NLSE",
    progress_message = "NLSE",
    kwargs...
    )

    progress ? (prog = Juno.ProgressBar(name=progress_name)) : prog = nothing

    cache = alg_cache(alg, u, rate_prototype, uEltypeNoUnits, tTypeNoUnits,
                      uprev, uprev2, f, t, dt, reltol_internal, Val{isinplace(prob)})

    sol = build_solution(prob, alg, ts, timeseries,
                         dense=dense, k=ks,
                         interp=id, calculate_error = false)

    u = copy(prob.u0)
    uType = typeof(u)
    uprev = copy(u)
    t = prob.tspan[1]
    tType = typeof(t)
    CacheType = typeof(cache)

    integrator = NLSEIntegrator{algType, uType, tType, CacheType}(
        sol, u, uprev, t, dt, alg, cache)

    initialize!(integrator, integrator.cache)

    integrator
end

function solve!(integrator::NLSEIntegrator)
    @inbounds while !isempty(integrator.opts.tstops)
        integrator.t = pop!(integrator.opts.tstops)
        perform_step!(integrator, integrator.cache)
        push!(integrator.out, integrator.u)
    end

    build_solution(p, integrator.alg, integrator.tstops, integrator.out)

    if typeof(integrator.sol.prob.f) <: Tuple
        f = integrator.sol.prob.f[1]
    else
        f = integrator.sol.prob.f
    end

    if has_analytic(f)
        calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
    end
    integrator.sol.retcode = :Success
    nothing
end
