# Simulate suprcontinuum generation for parameters similar to Fig. 3 of
# Dudley et al, RMP 78 1135 (2006)
# Parameters from www.scgbook.info
T_T0 = 0.0284ps
T_FWHM = T_T0 * 2*log(1+sqrt(2)) # pulse time paramater passed as FWHM
betacoeff = [0.0/m, 0.0ps/m, -11.83e-3ps^2/m, 8.1038e-5ps^3/m, -9.5205e-8ps^4/m,
             2.0737e-10ps^5/m, -5.3943e-13ps^6/m,
             1.13486e-15ps^7/m, -2.5495e-18ps^8/m, 3.0524e-21ps^9/m, -1.7140e-24ps^10/m]
γnl = 0.11W/m
ω0 = 2π*c0/(835nm) |> THz# c is m/ps, exported from NonlinearOptics
P_pulse = 100.0W
fiber_length = 0.15m
dudley = ToyModel(;ω0=ω0, betacoeff=betacoeff,
                   nonlinearcoeff=γnl,
                   linearloss=0.0/m,
                   power_in=P_pulse,
                   pulsetime=T_FWHM,
                   length=fiber_length,
                   has_shock=true,
                   has_raman=true)

sol_dudley = solve(dudley, DynamicNLSE(); time_pts=2^13, time_window=12.5ps)
# solves in 1.3 seconds, compared to 27 seconds for scgbook code!
