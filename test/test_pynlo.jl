T_FWHM = 0.05ps # pulse time paramater passed as FWHM
betacoeff = [0.0/m, 0.0ps/m, -0.12ps^2/m, 0.0ps^3/m, 0.005e-3ps^4/m]
γnl = 1.0W/m
ω0 = 2π*c0/(1550nm) |> THz
P_pulse = 10.0W
fiber_length = 0.020m
nlo = ToyModel(;ω0=ω0, betacoeff=betacoeff,
                   nonlinearcoeff=γnl,
                   linearloss=0.0/m,
                   power_in=P_pulse,
                   pulsetime=T_FWHM,
                   length=fiber_length,
                   has_shock=true,
                   has_raman=true)

sol_nlo = solve(nlo, DynamicNLSE();
                tpoints=2^13, time_window=10.0ps)
