function solve{algType<:AbstractNLSEAlgorithm}(
    prob::AbstractNLSEProblem,
    alg::algType,
    timeseries=[],
    ts=[],
    kwargs...)

    integrator = init(prob, alg, timeseries, ts; kwargs...)
    solve!(integrator)
    integrator.sol
end

function init{algType<:AbstractNLSEAlgorithm}(
    prob::AbstractNLSEAlgorithm,
    alg::algType;
    timeseries_init=typeof(prob.u0)[],
    timeseries_steps = 1,
    saveat = eltype(prob.tspan)[],
    tstops = eltype(prob.tspan)[],
    save_everystep = isempty(saveat),
    save_timeseries = nothing,
    save_start = true,
    dense = save_everystep && !(typeof(alg) <: Discrete),
    dt = typeofalg <: Discrete && isempty(tstops) ? eltype(prob.tspan)(1) : eltype(prob.tspan)(0),
    adaptive = isadaptive(alg),
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
