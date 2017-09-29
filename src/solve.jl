function solve(
    prob::NLSEProblem,
    alg::NLSESplitStepFourier;
    solver::Symbol=:Direct,
    save_everystep::Bool=false,
    timeseries_steps::Int = 100,
    autodiff::Bool=false,
    method=:trust_region,
    show_trace=false,
    iterations=1000,
    progress_steps::Int=1000,
    progressbar::Bool=false,
    progressbar_name="NLSE",
    kwargs...
    )

    #Unroll important constants
    @unpack f, u0, Du, analytic, numvars = prob

    if dt != 0
        numiters = round(Int64, tspan[2]/dt)
    else
        numiters = 0
    end

    #Set initial
    u = copy(u0)
    t = tspan[1]

    #Setup timeseries
    timeseries = Vector{typeof(u)}(0)
    push!(timeseries, copy(u))
    ts = Float64[t]

    sqrtdt = sqrt(dt)

    #Split-Step loop
    u, timeseries, ts = nlse_solve(NLSEIntegrator)
end
