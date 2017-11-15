function build_problem(model::Model, probtype;
                       tpoints=2^13, time_window=10.0, kwargs...)

    ω = get_ωmesh(tmesh)
    const dt_mesh = tmesh[2]-tmesh[1]
    const ω = fftshift(ω)
    const ω0 = model.ω0

    const dispersion = get_func(model.dispersion)
    const nonlinearcoeff = get_func(model.nonlinearcoeff)
    const α = get_func(model.linearloss)

    const D = @. 1im * dispersion(ω) - α(ω)/2
    const γnl = @. nonlinearcoeff(ω)

    const has_raman = model.has_raman
    const has_shock = model.has_shock

    const shock_term = build_model_shock(model, ω)
    const frac_raman = 0.18
    const raman_response = build_model_raman(model, tmesh)


    function f(z, u, du)
        @. du = u * exp(D * z)
        planned_fft! * du
        if has_raman
            raman_term = planned_ifft! * copy(abs2.(du)) # prepare intensity for convolution
            # mulitiply in fourier space for convolution
            @. raman_term = raman_term * raman_response
            planned_fft! * raman_term # fourier transform back
            @. du = du * ((1-frac_raman)*abs2(du) + frac_raman*raman_term)
        else
            @. du = du * abs2(du)
        end
        planned_ifft! * du
        @. du = 1im * γnl * shock_term * du * exp(-D * z)
    end

    if typeof(probtype) <: DynamicIkeda
        ikeda_callback = build_ikedacallback(model)
        prob = ODEProblem(f,
                          planned_ifft! * u0,
                          zspan;
                          kwargs...)
        prob_exp = DynamicIkedaProblem(prob,
                                       fftshift(ω),
                                       ω0,
                                       tmesh,
                                       planned_fft!,
                                       planned_ifft!,
                                       D,
                                       ikeda_callback)
    else
        prob = ODEProblem(f,
                          planned_ifft! * u0,
                          zspan;
                          kwargs...)
        prob_exp = DynamicNLSEProblem(prob,
                                  fftshift(ω),
                                  ω0,
                                  tmesh,
                                  planned_fft!,
                                  planned_ifft!,
                                  D)
    end
    return prob_exp
end

# takes a scalar and turns it into a callable
function get_func(attr::Number)
    f(x...) = attr
    return f
end
get_func(attr) = attr

function build_problem(model::ToyModel, probtype::Union{DynamicNLSE, DynamicIkeda};
                   tpoints=2^10, time_window=10.0, kwargs...)
    const α = model.linearloss
    const γnl = model.nonlinearcoeff
    const L = model.length
    beta = model.betacoeff
    const beta_coeff = beta .* [1/factorial(n) for n=0:(length(beta)-1)]
    tmesh = linspace(-time_window/2, time_window/2, tpoints)
    zspan = (0.0, model.length)

    ω = get_ωmesh(tmesh)
    const dt_mesh = tmesh[2]-tmesh[1]
    const ω = fftshift(ω)
    const ω0 = model.ω0

    u0 = derive_pulse(model.power_in, model.pulsetime).(tmesh)
    const planned_fft! = plan_fft!(u0, flags=FFTW.MEASURE)
    const planned_ifft! = plan_ifft!(u0, flags=FFTW.MEASURE)
    const D = 1im * Poly(beta_coeff, :ω).(ω) - α/2

    const has_raman = model.has_raman
    const has_shock = model.has_shock

    const shock_term = build_model_shock(model, ω)
    const frac_raman = 0.18
    const raman_response = build_model_raman(model, tmesh)

    function f(z, u, du)
        @. du = u * exp(D * z)
        planned_fft! * du
        if has_raman
            raman_term = planned_ifft! * copy(abs2.(du)) # prepare intensity for convolution
            # mulitiply in fourier space for convolution
            @. raman_term = raman_term * raman_response
            planned_fft! * raman_term # fourier transform back
            @. du = du * ((1-frac_raman)*abs2(du) + frac_raman*raman_term)
        else
            @. du = du * abs2(du)
        end
        planned_ifft! * du
        @. du = 1im * γnl * shock_term * du * exp(-D * z)
    end

    if typeof(probtype) <: DynamicIkeda
        ikeda_callback = build_ikedacallback(model)
        prob = ODEProblem(f,
                          planned_ifft! * u0,
                          zspan;
                          kwargs...)
        prob_exp = DynamicIkedaProblem(prob,
                                       fftshift(ω),
                                       ω0,
                                       tmesh,
                                       planned_fft!,
                                       planned_ifft!,
                                       D,
                                       ikeda_callback)
    else
        prob = ODEProblem(f,
                          planned_ifft! * u0,
                          zspan;
                          kwargs...)
        prob_exp = DynamicNLSEProblem(prob,
                                  fftshift(ω),
                                  ω0,
                                  tmesh,
                                  planned_fft!,
                                  planned_ifft!,
                                  D)
    end
    prob_exp
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

    const has_raman = model.has_raman
    const has_shock = model.has_shock

    const shock_term = build_model_shock(model, ω)
    const frac_raman = 0.18
    const raman_response = build_model_raman(model, tmesh)

    const planned_fft! = plan_fft!(u0, flags=FFTW.MEASURE)
    const planned_ifft! = plan_ifft!(u0, flags=FFTW.MEASURE)
    const D = FSR * (-1im * L * Poly(beta_coeff, :ω).(ω) - α - 1im * detuning)
    #TODO: Macro for adding terms to function
    function f(z, u, du)
        @. du = (u * exp(D * z))
        planned_fft! * du
        if has_raman
            raman_term = planned_ifft! * copy(abs2.(du)) # prepare intensity for convolution
            # mulitiply in fourier space for convolution
            @. raman_term = raman_term * raman_response
            planned_fft! * raman_term # fourier transform back
            @. du = du * ((1-frac_raman)*abs2(du) + frac_raman*raman_term)
        else
            @. du = du * abs2(du)
        end
        planned_ifft! * du
        @. du = (1im * γnl * L * FSR * du) * exp(-D * z)
        du[1] = sqrtcoupling * Ein * FSR * exp(-D[1] * z)
    end

    prob = ODEProblem(f, planned_ifft! * u0, τspan; kwargs...)
    prob_LL = DynamicLLProblem(prob, fftshift(ω), ω0, tmesh, planned_fft!, planned_ifft!, D)
end

function build_ikedacallback(model)
    L = model.length
    FSR = model.FSR
    const Ein = sqrt(model.power_in)
    const coupling = sqrt(model.coupling)
    const transmission = sqrt(1-coupling^2)
    const phase_shift = 2pi*L-model.detuning
    # summation is in fourier space
    function affect!(integrator)
        integrator.u = transmission * exp(1im*phase_shift) * integrator.u
        integrator.u[1] += coupling * Ein
    end

    ikeda_callback = PeriodicCallback(affect!, L)
end

function build_model_shock(model, ω)
    ω0 = model.ω0

    if model.has_shock
        shock_term = (ω + ω0)/ω0
    else
        shock_term = one(ω0)
    end
    shock_term
end

function build_model_raman(model, tmesh)
    raman_response = (1+0im)*similar(tmesh)
    if model.has_raman
        tau1 = 0.0122; tau2 = 0.032
        raman_timeresponse = @. (1+0im)*(tmesh > 0) * (tau1^2 + tau2^2)/tau1/tau2^2*exp(-tmesh/tau2)*sin(tmesh/tau1)
        raman_response = length(tmesh)*ifft(fftshift(raman_timeresponse))
    end
    raman_response
end
