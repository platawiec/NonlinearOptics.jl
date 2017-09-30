function solve(
    prob::AbstractNLSEProblem{uType, tType, zType},
    alg::NLSESplitStepFourier,
    timeseries=[],
    ts=[],
    kwargs...)

    integrator = init(prob, alg, timeseries, ts; kwargs...)
    solve!(integrator)
    integrator.sol
end

function init(
    prob::AbstractNLSEProblem{uType, tType, zType},
    alg::NLSESplitStepFourier,
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
    initialize_integrator = true,
    kwargs...
    )

    tType = eltype(prob.tspan)
    tspan = prob.tspan
    tdir = sign(tspan[end]-tspan[1])

    t = tspan[1]

    if dt == tType(0) && isempty(tstops)
        error("NLSE requires a choice of dt or choosing the tstops")
    end

    if tspan[1] == tspan[end]
        error("Timespan is trivial")
    end

    u = recursivecopy(prob.u0)
    uType = typeof(u)
    uElType = recursive_eltype(u)


    integrator = NLSEIntegrator()

    if initialize_integrator
        initialize!(integrator)
        #initialize!(callbacks_internal, t, u, integrator)
    end

    integrator
end
