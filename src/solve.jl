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
    adaptive = false,
    abstol = nothing,
    reltol = nothing,
    maxiters = 1000000,
    verbose = true,
    progress = false,
    progress_steps = 1000,
    progress_name = "NLSE",
    progress_message = "NLSE",
    )
end

function solve!(integrator::NLSEIntegrator)
    @inbounds while !isempty(integrator.opts.tstops)
        while integrator.tdir*integrator.t < integrator.tdir*top(integrator.opts.tstops)
            loopheader!(integrator)
            perform_step!(integrator,integrator.cache)
            loopfooter!(integrator)
            if isempty(integrator.opts.tstops)
                break
            end
        end
        handle_tstop!(integrator)
    end
    postamble!(integrator)

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
