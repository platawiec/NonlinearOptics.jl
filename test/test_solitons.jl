# Second-order Soliton model from J. Lightwave Tech. Vol. 25 No. 12 Dec. 2007
# Johan Hult
T_FWHM = 0.1ps
T_T0 = T_FWHM/(2*log(1+sqrt(2)))
β₂ = -0.01ps^2/m
γnl = 0.01W/m
P_pulse = (abs(β₂/γnl)/(T_T0))^2
soliton_period = π*T_T0^2/(2*abs(β₂))
soliton_fundamental = ToyModel(;betacoeff=[0.0/m, 0.0ps/m, β₂],
                                nonlinearcoeff=γnl,
                                linearloss=0.0/m,
                                power_in=P_pulse,
                                pulsetime=T_FWHM,
                                length=soliton_period)

P_pulse_secondorder = (abs(β₂/γnl)/(T_T0)*2)^2
soliton_secondorder = ToyModel(;ω0=200.0THz, betacoeff=[0.0/m, 0.0ps/m, β₂],
                                nonlinearcoeff=γnl,
                                linearloss=0.0/m,
                                power_in=P_pulse_secondorder,
                                pulsetime=T_FWHM,
                                length=soliton_period)

sol_fundamental = solve(soliton_fundamental, DynamicNLSE();
                        tpoints=2^12, time_window=25.0ps)
sol_secondorder = solve(soliton_secondorder, DynamicNLSE();
                        tpoints=2^12, time_window=25.0ps)

sq_diff(x, y) = abs(sum(abs2.(x) - abs2.(y)))
@test sq_diff(sol_fundamental(0m), sol_fundamental(soliton_period)) < 15
@test sq_diff(sol_secondorder(soliton_period), sol_secondorder(0m)) < 15
@test sq_diff(sol_fundamental(0m), sol_secondorder(0m)) > 5
@test sq_diff(FT(sol_fundamental, 0m), FT(sol_fundamental, soliton_period)) < 15/ps^2
@test sq_diff(sol_fundamental(soliton_period/2), sol_fundamental(0m)) < 15
