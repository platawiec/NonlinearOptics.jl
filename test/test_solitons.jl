# Second-order Soliton model from J. Lightwave Tech. Vol. 25 No. 12 Dec. 2007
# Johan Hult
T_FWHM = 0.1
T_T0 = T_FWHM/(2*log(1+sqrt(2)))
β₂ = -0.01
γnl = 0.01
P_pulse = (abs(β₂/γnl)/(T_T0))^2
soliton_period = π*T_T0^2/(2*abs(β₂))
soliton_fundamental = ToyModel(;betacoeff=[0.0, 0.0, β₂],
                                nonlinearcoeff=γnl,
                                linearloss=0.0,
                                power_in=P_pulse,
                                pulsetime=T_FWHM,
                                length=soliton_period)

P_pulse_secondorder = (abs(β₂/γnl)/(T_T0)*2)^2
soliton_secondorder = ToyModel(;betacoeff=[0.0, 0.0, β₂],
                                nonlinearcoeff=γnl,
                                linearloss=0.0,
                                power_in=P_pulse_secondorder,
                                pulsetime=T_FWHM,
                                length=soliton_period)

prob_fundamental = build_problem(soliton_fundamental, DynamicNLSE();
                                tpoints=2^12, time_window=25.0)
prob_secondorder = build_problem(soliton_secondorder, DynamicNLSE();
                                tpoints=2^12, time_window=25.0)

sol_fundamental = solve(prob_fundamental)
sol_secondorder = solve(prob_secondorder)

sq_diff(x, y) = abs(sum(abs2.(x) - abs2.(y)))
@test sq_diff(sol_fundamental(0), sol_fundamental(soliton_period)) < 5
@test sq_diff(sol_secondorder(soliton_period), sol_secondorder(0)) < 5
@test sq_diff(sol_fundamental(0), sol_secondorder(0)) > 5
@test sq_diff(FT(sol_fundamental, 0), FT(sol_fundamental, soliton_period)) < 5
@test sq_diff(sol_fundamental(soliton_period/2), sol_fundamental(0)) < 5
