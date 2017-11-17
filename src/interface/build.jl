function build_problem(model, probtype::DynamicNLSE;
                       tpoints=2^13, time_window=10.0, kwargs...)

    tmesh = linspace(-time_window/2, time_window/2, tpoints)
    zspan = (0.0, pathlength(model))

    ω = get_ωmesh(tmesh)
    const dt_mesh = tmesh[2]-tmesh[1]
    const ω = fftshift(ω)
    const ω0 = get_ω0(model)
    u0 = derive_pulse.(model, tmesh)

    f, D, planned_fft!, planned_ifft! = build_GNLSE(model, ω, tmesh)
    prob = ODEProblem(
        f,
        planned_ifft! * u0,
        zspan;
        kwargs...
    )
    prob_exp = DynamicNLSEProblem(
        prob,
        fftshift(ω),
        ω0,
        tmesh,
        planned_fft!,
        planned_ifft!,
        D
    )
    return prob_exp
end

function build_problem(model, probtype::DynamicIkeda;
                       tpoints=2^13, time_window=10.0, kwargs...)

    tmesh = linspace(-time_window/2, time_window/2, tpoints)
    zspan = (0.0, pathlength(model))

    ω = get_ωmesh(tmesh)
    const dt_mesh = tmesh[2]-tmesh[1]
    const ω = fftshift(ω)
    const ω0 = model.ω0
    u0 = derive_pulse.(model, tmesh)

    f, D, planned_fft!, planned_ifft! = build_GNLSE(model, ω, tmesh)
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
    return prob_exp
end

# takes a scalar and turns it into a callable
function get_func(attr::Number)
    f(x...) = attr
    return f
end
get_func(attr) = attr

function build_problem(model::ToyModel, ::DynamicLL;
                       tpoints=2^10, time_window=10.0, kwargs...)
    FSR = model.FSR
    tmesh = linspace(-1/2/FSR, 1/2/FSR, tpoints)
    τspan = (0.0, time_window)

    ω = get_ωmesh(tmesh)
    ω = fftshift(ω)
    ω0 = 0.0#toy model ω0 is 0

    #TODO: random start conditions don't seem to do anything
    u0 = derive_pulse(model.power_in, model.pulsetime).(tmesh)

    f, D, planned_fft!, planned_ifft! = build_GLLE(model, ω, tmesh)

    prob = ODEProblem(f, planned_ifft! * u0, τspan; kwargs...)
    prob_LL = DynamicLLProblem(prob, fftshift(ω), ω0, tmesh, planned_fft!, planned_ifft!, D)
end

function build_GNLSE(model::ToyModel, ω, tmesh)

    const α = linearloss(model, ω)
    const γnl = nonlinearcoeff(model, ω)
    const beta = stationarydispersion(model, ω)

    const planned_fft! = plan_fft!(Vector((1+0im)*tmesh), flags=FFTW.MEASURE)
    const planned_ifft! = plan_ifft!(Vector((1+0im)*tmesh), flags=FFTW.MEASURE)
    const D = 1im .* beta .- α/2

    const use_raman = has_raman(model)

    const shock_term = build_model_shock(model, ω)
    const raman_response = build_model_raman(model, tmesh)

    const frac_raman = 0.18

    function f(z, u, du)
        @. du = u * exp(D * z)
        planned_fft! * du
        if use_raman
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
    return f, D, planned_fft!, planned_ifft!
end

function build_GNLSE(model::Model, ω, tmesh)

    const α = linearloss(model, ω)
    const γnl = nonlinearcoeff(model, ω)
    const beta = stationarydispersion(model, ω)

    const planned_fft! = plan_fft!(Vector((1+0im)*tmesh), flags=FFTW.MEASURE)
    const planned_ifft! = plan_ifft!(Vector((1+0im)*tmesh), flags=FFTW.MEASURE)
    const D = 1im .* beta .- α/2

    const use_raman = has_raman(model)

    const shock_term = build_model_shock(model, ω)
    const raman_response = build_model_raman(model, tmesh)

    function f(z, u, du)
        @. du = u * exp(D * z)
        planned_fft! * du
        if use_raman
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
    return f, D, planned_fft!, planned_ifft!
end

function build_GLLE(model, ω, tmesh)
    FSR = model.FSR
    α = linearloss(model, ω)
    γnl = nonlinearcoeff(model, ω)
    L = model.length
    detuning = model.detuning
    sqrtcoupling = sqrt(model.coupling)
    Ein = sqrt(model.power_in)
    beta = model.betacoeff
    beta_coeff = beta .* [(1im)^n/factorial(n) for n=0:(length(beta)-1)]
    const has_raman = model.has_raman
    const has_shock = model.has_shock

    const shock_term = build_model_shock(model, ω)
    const frac_raman = 0.18
    const raman_response = build_model_raman(model, tmesh)

    const planned_fft! = plan_fft!(Vector((1+0im)*tmesh), flags=FFTW.MEASURE)
    const planned_ifft! = plan_ifft!(Vector((1+0im)*tmesh), flags=FFTW.MEASURE)
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
    return f, D, planned_fft!, planned_ifft!
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

function build_model_shock(model::ToyModel, ω)
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
