T_FWHM = 0.05 # pulse time paramater passed as FWHM
betacoeff = [0., 0., -120.0, 0.0, 0.005]
γnl = 1.
ω0 = 2π*c/(1550e-9) # c is m/ps, exported from NonlinearOptics
P_pulse = 10000.
fiber_length = 0.020
nlo = ToyModel(;ω0=ω0, betacoeff=betacoeff,
                   nonlinearcoeff=γnl,
                   linearloss=0.0,
                   power_in=P_pulse,
                   pulsetime=T_FWHM,
                   length=fiber_length,
                   has_shock=false,
                   has_raman=false)

prob_nlo = build_problem(nlo, DynamicNLSE();
                            tpoints=2^13, time_window=10.)
sol_nlo = solve(prob_nlo)
