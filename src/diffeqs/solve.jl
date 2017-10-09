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
    ts_init=eltype(prob.tspan)[], zmesh_init=eltype(prob.zmesh)[],
    ks_init=[];
    timeseries_steps = 1,
    saveat = eltype(prob.tspan)[],
    tstops = eltype(prob.tspan)[],
    save_everystep = isempty(saveat),
    save_timeseries = nothing,
    save_start = true,
    dense = save_everystep,
    dt = (prob.tspan[2]-prob.tspan[1])/100,
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
    zs = prob.zmesh
    dktilde = 2pi/maximum(prob.zmesh)
    ktilde_max = dktilde/2.0*(length(prob.zmesh)-1)
    ktilde = collect(-ktilde_max:dktilde:ktilde_max)
    ktilde = fftshift(ktilde)

    ks = ks_init
    timeseries = convert(Vector{uType}, timeseries_init)
    ts = convert(Vector{tType}, ts_init)

    cache = alg_cache(alg, u)

    CacheType = typeof(cache)

    sol = build_solution(prob, alg, ts, zs, timeseries)

    integrator = NLSEIntegrator{algType, uType, tType, typeof(tstops), typeof(ktilde), CacheType, typeof(prob.N), typeof(prob.D)}(
        prob.N, prob.D, sol, u, uprev, t, dt, tstops, ktilde, alg, cache)

    initialize!(integrator, integrator.cache)

    return integrator
end

function solve!(integrator::NLSEIntegrator)
    @inbounds while !isempty(integrator.tstops)
        integrator.t = shift!(integrator.tstops)
        perform_step!(integrator, integrator.cache)
        push!(integrator.sol.u, integrator.u)
        push!(integrator.sol.t, integrator.t)
    end

    integrator.sol.retcode = :Success
    nothing
end
