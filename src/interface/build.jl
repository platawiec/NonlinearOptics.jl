function build_problem(model::ToyModel, ::DynamicNLSE;
                       tpoints=2^10, time_window=10.0, kwargs...)
    α = model.linearloss
    γnl = model.nonlinearcoeff
    L = model.length
    beta = model.betacoeff
    beta_coeff = beta .* [(1im)^n/factorial(n) for n=0:(length(beta)-1)]
    tmesh = linspace(-time_window/2, time_window/2, tpoints)
    zspan = (0.0, model.length)

    dt_mesh = tmesh[2]-tmesh[1]
    ω = collect(2pi*(-length(tmesh)/2:length(tmesh)/2-1)/(length(tmesh)*dt_mesh))
    ω = fftshift(ω)
    ω0 = 0.0#toy model ω0 is 0

    u0 = derive_pulse(model.power_in, model.pulsetime).(tmesh)
    planned_fft = plan_fft(u0, flags=FFTW.MEASURE)
    planned_ifft = plan_ifft(u0, flags=FFTW.MEASURE)
    D = -1im * Poly(beta_coeff, :ω).(ω) - α/2
    #TODO: Macro for adding terms to function
    function f(z, u)
        uT = planned_fft * (u .* exp(D * z))
        1im*γnl.*planned_ifft * (uT.*abs2.(uT)) .* exp(-D * z)
    end

    prob = ODEProblem(f, planned_ifft * u0, zspan; kwargs...)
    prob_NLSE = DynamicNLSEProblem(prob, fftshift(ω), ω0, tmesh, planned_fft, planned_ifft, D)
end
