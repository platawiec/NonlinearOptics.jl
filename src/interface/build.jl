function build_problem(model::ToyModel, ::DynamicNLSE;
                       tpoints=2^10, time_window=10.0, kwargs...)
    α = model.linearloss
    γnl = model.nonlinearcoeff
    L = model.length
    Ein = sqrt(model.power_in)
    detuning = model.detuning
    sqrtcoupling = sqrt(model.coupling)
    beta = model.betacoeff
    beta_coeff = beta .* [(1im)^n/factorial(n) for n=0:(length(beta)-1)]

    tmesh = linspace(-model.pulsetime*10, model.pulsetime*10, tpoints)
    dt_mesh = tmesh[2]-tmesh[1]
    const sqrtdt = sqrt(dt_mesh)
    zspan = (0.0, model.length)

    dω = 2pi/maximum(tmesh)
    ω_max = dω/2.0*(length(tmesh)-1)
    ω = collect(-ω_max:dω:ω_max)
    ω = fftshift(ω)
    ω0 = 0.0#toy model ω0 is 0

    u0 = derive_pulse(Ein, 1.0, model.pulsetime/(2*log(1+sqrt(2))))(tmesh)
    planned_fft = plan_fft(u0)
    planned_ifft = plan_ifft(u0)
    D = 1im * Poly(beta_coeff, :ω).(ω) - α/2
    #TODO: Macro for adding terms to function
    function f(z, u)
        uT = planned_fft * (u .* exp(D * z))
        1im*γnl.*planned_ifft * (uT.*abs2.(uT)) .* exp(-D .* z)
    end

    prob = ODEProblem(f, planned_ifft * u0, zspan)
    prob_NLSE = DynamicNLSEProblem(prob, fftshift(ω), ω0, tmesh, planned_fft, planned_ifft, D)
end
