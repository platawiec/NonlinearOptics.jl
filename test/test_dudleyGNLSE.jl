# Simulate suprcontinuum generation for parameters similar to Fig. 3 of
# Dudley et al, RMP 78 1135 (2006)
# Parameters from www.scgbook.info
T_T0 = 0.0284
T_FWHM = T_T0 * 2*log(1+sqrt(2)) # pulse time paramater passed as FWHM
betacoeff = [0., 0., -11.83e-3, 8.1038e-5, -9.5205e-8, 2.0737e-10, -5.3943e-13,
             1.13486e-15, -2.5495e-18, 3.0524e-21, -1.7140e-24]
γnl = 0.11
ω0 = 2π*c/(835e-9) # c is m/ps, exported from NonlinearOptics
P_pulse = 100.
fiber_length = 0.15
dudley = ToyModel(;ω0=ω0, betacoeff=betacoeff,
                   nonlinearcoeff=γnl,
                   linearloss=0.0,
                   power_in=P_pulse,
                   pulsetime=T_FWHM,
                   length=fiber_length,
                   has_shock=true,
                   has_raman=true)

prob_dudley = build_problem(dudley, DynamicNLSE();
                            tpoints=2^13, time_window=12.5)
sol_dudley = solve(prob_dudley)
# solves in 1.3 seconds, compared to 27 seconds for scgbook code!
