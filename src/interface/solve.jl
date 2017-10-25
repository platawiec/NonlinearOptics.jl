function solve(model::Model, exp::Union{DynamicLL, DynamicNLSE}; kwargs...)
    prob = build_problem(model.laser, model.structure, exp; kwargs...)

    solve(prob, SymmetrizedSplitStep(); kwargs...)
end

function solve(model::Model, exp::DynamicIkeda;
               steps_per_roundtrip=100, kwargs...)
    prob = build_problem(model.laser, model.structure, exp; kwargs...)
    ikeda_callback = build_ikedacallback(model.laser, model.structure)

    dz = circumference(model.structure)/steps_per_roundtrip
    solve(prob, SymmetrizedSplitStep(), callback=ikeda_callback, dt=dz; kwargs...)
end

function solve(model::ToyModel, exp::DynamicLL; tpoints=2^10, time_window=100, kwargs...)
    FSR = model.FSR
    α = model.linearloss
    γnl = model.nonlinearcoeff
    L = model.length
    Ein = sqrt(model.power_in)
    detuning = model.detuning
    sqrtcoupling = sqrt(model.coupling)
    beta = model.betacoeff
    beta_coeff = beta .* [(1im)^n/factorial(n) for n=0:(length(beta)-1)]

    tmesh = linspace(-1/2/FSR, 1/2/FSR, tpoints)
    dt_mesh = tmesh[2]-tmesh[1]
    const sqrtdt = sqrt(dt_mesh)
    τspan = (0.0, time_window)

    u0 = derive_pulse(Ein, 1.0, 1/FSR*0.01)
    N(z, u) = (1im*γnl*L*FSR*abs2(u) - FSR*(α+1im*detuning))*u + FSR*sqrtcoupling*Ein
    D(ω, u) = -1im * FSR * L * Poly(beta_coeff, :ω)(ω)

    prob = NLSEProblem(N, D, u0, τspan, tmesh)
    solve(prob, SymmetrizedSplitStep(); kwargs...)
end

"""
    build_problem()
Sets up a Lugatio-Lefever problem
"""
function build_problem(laser::CWLaser, res::AbstractResonator, ::DynamicLL;
                       dispersion_order=2, additional_terms=[],
                       time_window=0.01, tpoints=2^10, kwargs...)
    # TODO: support mode coupling
    # TODO: currently only supports one mode
    mode = res.modes[1]

    FSR = get_FSR(res, mode, laser.frequency)
    α = mode.linearloss(laser.frequency)
    γnl = get_nonlinearcoeff(res, mode, laser.frequency)
    L = circumference(res)
    Ein = sqrt(laser.power)
    detuning = laser.detuning
    const sqrtcoupling = sqrt(mode.coupling(laser.frequency))
    beta = get_beta(mode, laser.frequency, dispersion_order)
    beta_coeff = beta .* [(1im)^n/factorial(n) for n=0:(length(beta)-1)]
    beta_coeff[1] = 0
    beta_coeff[2] = 0

    tmesh = linspace(-1/2/FSR, 1/2/FSR, tpoints)
    dt_mesh = tmesh[2]-tmesh[1]
    const sqrtdt = sqrt(dt_mesh)
    τspan = (0.0, time_window)

    # random definition
    # u0(t) = ((1+0im)*rand()+(0+1im)*rand())/sqrtdt
    u0 = derive_pulse(Ein, 1.0, 1/FSR*0.1)
    if :self_steepening in additional_terms
        ω0 = getω(laser.frequency)
        self_steepening = :(- γnl/ω0 * diff_cyclic(abs2(u)) / dt_mesh * u)
    end
    if :raman_response in additional_terms
        raman_response = :(1im*γnl*material_raman_response(wg.material))
    end

    #TODO: Macro for adding terms to function
    N(z, u) = (1im*γnl*L*FSR*abs2(u) - FSR*(α+1im*detuning))*u + FSR*sqrtcoupling*Ein
    D(ω, u) = -1im * FSR * L * Poly(beta_coeff, :ω)(ω)

    prob = NLSEProblem(N, D, u0, τspan, tmesh)
end

function build_problem(laser::PulsedLaser, wg::Waveguide, ::DynamicNLSE;
                       dispersion_order=2, additional_terms=[],
                       time_window=2.0, tpoints=2^10, kwargs...)
    # TODO: support mode coupling
    # TODO: currently only supports one mode
    mode = wg.modes[1]

    α = mode.linearloss(laser)
    γnl = get_nonlinearcoeff(wg, mode, laser.frequency)
    L = wg.length

    beta = get_beta(mode, laser.frequency, dispersion_order)
    beta_coeff = beta .* [(1im)^n/factorial(n) for n=0:(length(beta)-1)]
    beta_coeff[1] = 0
    beta_coeff[2] = 0
    tmesh = linspace(-time_window/2, time_window/2, tpoints)
    dt_mesh = tmesh[2] - tmesh[1]
    zspan = (0.0, L)
    const coupling = mode.coupling(laser)

    u0(t) = coupling * laser.pulse_init(t)
    if :self_steepening in additional_terms
        ω0 = getω(laser.frequency)
        self_steepening = :(- γnl/ω0 * diff_cyclic(abs2(u)) / dt_mesh * u)
    end
    if :raman_response in additional_terms
        raman_response = :(1im*γnl*material_raman_response(wg.material))
    end

    #TODO: Macro for adding terms to function
    N(z, u) = (-α/2 + 1im * γnl * abs2(u)) * u
    D(ω, u) = -1im * Poly(beta_coeff, :ω)(ω)

    prob = NLSEProblem(N, D, u0, zspan, tmesh)
end

function build_problem(laser::CWLaser, res::AbstractResonator, ::DynamicIkeda;
                       dispersion_order=2, additional_terms=[],
                       round_trips=100, tpoints=2^10, kwargs...)
    # TODO: support mode coupling
    # TODO: currently only supports one mode
    mode = res.modes[1]

    α = mode.linearloss(laser)
    γnl = get_nonlinearcoeff(res, mode, laser.frequency)
    L = circumference(res)

    FSR = get_FSR(res, mode, laser.frequency)
    beta = get_beta(mode, laser.frequency, dispersion_order)
    beta_coeff = beta .* [(1im)^n/factorial(n) for n=0:(length(beta)-1)]
    beta_coeff[1] = 0
    beta_coeff[2] = 0
    tmesh = linspace(-1/2/FSR, 1/2/FSR, tpoints)
    dt_mesh = tmesh[2] - tmesh[1]
    zspan = (0.0, L*round_trips)
    const Ein = sqrt(laser.power)
    const coupling = sqrt(mode.coupling(laser))
    const transmission = sqrt(1-coupling^2)
    const phase_shift = 2pi*L-laser.detuning

    u0 = derive_pulse(Ein, 1.0, 1/FSR*0.1)
    if :self_steepening in additional_terms
        ω0 = getω(laser.frequency)
        self_steepening = :(- γnl/ω0 * diff_cyclic(abs2(u)) / dt_mesh * u)
    end
    if :raman_response in additional_terms
        raman_response = :(1im*γnl*material_raman_response(res.material))
    end
    #TODO: Macro for adding terms to function
    N(z, u) = (-α/2 + 1im * γnl * abs2(u)) * u
    D(ω, u) = -1im * Poly(beta_coeff, :ω)(ω)

    prob = NLSEProblem(N, D, u0, zspan, tmesh)
    prob
end

function build_ikedacallback(laser::CWLaser, res::AbstractResonator)
    mode = res.modes[1]

    L = circumference(res)
    FSR = get_FSR(res, mode, laser.frequency)
    const Ein = sqrt(laser.power)
    const coupling = sqrt(mode.coupling(laser))
    const transmission = sqrt(1-coupling^2)
    const phase_shift = 2pi*L-laser.detuning
    function affect!(integrator)
        integrator.u = transmission * exp(1im*phase_shift) * integrator.u + coupling * Ein
    end
    function condition(t, u, integrator)
        t % L
    end

    ikeda_callback = ContinuousCallback(condition, affect!)
end

function diff_cyclic(A)
    B = similar(A)
    @inbounds B[1] = A[1]-A[end]
    for i=2:length(B)
        @inbounds B[i] = A[i]-A[i-1]
    end
    B
end

function diff_cyclic!(A)
    @inbounds tmp = A[1]-A[end]
    for i=1:(length(A)-1)
        @inbounds A[i] = A[i+1]-A[i]
    end
    @inbounds A[end] = tmp
    A
end
"""
    build_problem()
Sets up a Lugatio-Lefever problem
"""
function build_problem(laser::CWLaser, res::AbstractResonator, ::DynamicLL;
                       dispersion_order=2, additional_terms=[],
                       time_window=0.01, tpoints=2^10, kwargs...)
    # TODO: support mode coupling
    # TODO: currently only supports one mode
    mode = res.modes[1]

    FSR = get_FSR(res, mode, laser.frequency)
    α = mode.linearloss(laser.frequency)
    γnl = get_nonlinearcoeff(res, mode, laser.frequency)
    L = circumference(res)
    Ein = sqrt(laser.power)
    detuning = laser.detuning
    const sqrtcoupling = sqrt(mode.coupling(laser.frequency))
    beta = get_beta(mode, laser.frequency, dispersion_order)
    beta_coeff = beta .* [(1im)^n/factorial(n) for n=0:(length(beta)-1)]
    beta_coeff[1] = 0
    beta_coeff[2] = 0

    tmesh = linspace(-1/2/FSR, 1/2/FSR, tpoints)
    dt_mesh = tmesh[2]-tmesh[1]
    const sqrtdt = sqrt(dt_mesh)
    τspan = (0.0, time_window)

    dω = 2pi/maximum(tmesh)
    ω_max = dω/2.0*(length(tmesh)-1)
    ω = collect(-ω_max:dω:ω_max)
    ω = fftshift(ω)
    # random definition
    # u0(t) = ((1+0im)*rand()+(0+1im)*rand())/sqrtdt
    u0 = derive_pulse(Ein, 1.0, 1/FSR*0.1)(tmesh)
    if :self_steepening in additional_terms
        ω0 = getω(laser.frequency)
        self_steepening = :(- γnl/ω0 * diff_cyclic(abs2(u)) / dt_mesh * u)
    end
    if :raman_response in additional_terms
        raman_response = :(1im*γnl*material_raman_response(wg.material))
    end

    planned_fft = plan_fft(u0)
    planned_ifft = plan_ifft(u0)
    #TODO: Macro for adding terms to function
    f(z, u) = (exp.(-z * planned_ifft * (-1im * FSR * L * Poly(beta_coeff, :ω).(ω) .* (planned_fft * u)))) .*
               (((1im*γnl*L*FSR*abs2.(u) - FSR*(α+1im*detuning)).*u + FSR*sqrtcoupling*Ein) .*
                exp.(z * planned_ifft * (-1im * FSR * L * Poly(beta_coeff, :ω).(ω) .* (planned_fft * u))))

    prob = ODEProblem(f, u0, τspan)
end

function solve(prob::AbstractNLOProblem, alg::AbstractODEAlgorithm=DP5(); kwargs...)
    sol = solve(prob.prob, alg; kwargs...)
    DynamicNLOSolution(sol, prob)
end
function solve(prob::DynamicIkedaProblem, alg::AbstractODEAlgorithm=DP5(); kwargs...)
    sol = solve(prob.prob, alg; callback=prob.ikeda_callback, kwargs...)
    DynamicNLOSolution(sol, prob)
end
