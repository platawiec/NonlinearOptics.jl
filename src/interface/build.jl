function build_problem(model::ToyModel, ::DynamicNLSE;
                       tpoints=2^10, time_window=10.0, kwargs...)
    α = model.linearloss
    γnl = model.nonlinearcoeff
    L = model.length
    beta = model.betacoeff
    beta_coeff = beta .* [(1im)^n/factorial(n) for n=0:(length(beta)-1)]
    tmesh = linspace(-time_window/2, time_window/2, tpoints)
    zspan = (0.0, model.length)

    ω = get_ωmesh(tmesh)
    ω = fftshift(ω)
    ω0 = 0.0#toy model ω0 is 0

    u0 = derive_pulse(model.power_in, model.pulsetime).(tmesh)
    planned_fft! = plan_fft!(u0, flags=FFTW.MEASURE)
    planned_ifft! = plan_ifft!(u0, flags=FFTW.MEASURE)
    D = -1im * Poly(beta_coeff, :ω).(ω) - α/2
    #TODO: Macro for adding terms to function
    function f(z, u, du)
        @. du = (u * exp(D * z))
        planned_fft! * du
        @. du = du * abs2(du)
        planned_ifft! * du
        @. du = 1im * γnl * du * exp(-D * z)
    end

    prob = ODEProblem(f, planned_ifft! * u0, zspan; kwargs...)
    prob_NLSE = DynamicNLSEProblem(prob, fftshift(ω), ω0, tmesh, planned_fft!, planned_ifft!, D)
end

function build_problem(model::ToyModel, ::DynamicLL;
                       tpoints=2^10, time_window=10.0, kwargs...)
    FSR = model.FSR
    α = model.linearloss
    γnl = model.nonlinearcoeff
    L = model.length
    detuning = model.detuning
    sqrtcoupling = sqrt(model.coupling)
    Ein = sqrt(model.power_in)
    beta = model.betacoeff
    beta_coeff = beta .* [(1im)^n/factorial(n) for n=0:(length(beta)-1)]
    tmesh = linspace(-1/2/FSR, 1/2/FSR, tpoints)
    τspan = (0.0, time_window)

    ω = get_ωmesh(tmesh)
    ω = fftshift(ω)
    ω0 = 0.0#toy model ω0 is 0

    if model.pulsetime == 0
        #TODO: random start conditions don't seem to do anything
        u0 = Ein + Ein*((1+0im)*rand(length(tmesh))+(0+1im)*rand(length(tmesh)))
    else
        u0 = derive_pulse(model.power_in, model.pulsetime).(tmesh)
    end
    planned_fft! = plan_fft!(u0, flags=FFTW.MEASURE)
    planned_ifft! = plan_ifft!(u0, flags=FFTW.MEASURE)
    D = FSR * (-1im * L * Poly(beta_coeff, :ω).(ω) - α - 1im * detuning)
    #TODO: Macro for adding terms to function
    function f(z, u, du)
        @. du = (u * exp(D * z))
        planned_fft! * du
        @. du = du * abs2(du) + sqrtcoupling * Ein / L / γnl /1im
        planned_ifft! * du
        @. du = (1im * γnl * L * FSR * du) * exp(-D * z)
    end

    prob = ODEProblem(f, planned_ifft! * u0, τspan; kwargs...)
    prob_LL = DynamicLLProblem(prob, fftshift(ω), ω0, tmesh, planned_fft!, planned_ifft!, D)
end
