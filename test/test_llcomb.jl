# Simple LL comb simulation
T_FWHM = 0.2
T_T0 = T_FWHM/(2*log(1+sqrt(2)))
β₂ = -0.013
γnl = 0.000032
α = 1.75e-5
coupling = 1.75e-5
P_in= 55.6e-3
FSR = 0.0182
L = 0.0119
detuning = 0.0012
comb_LL = ToyModel(;betacoeff=[0.0, 0.0, β₂],
                    nonlinearcoeff=γnl,
                    linearloss=α,
                    power_in=P_in,
                    pulsetime=T_FWHM,
                    length=L,
                    FSR=FSR,
                    detuning=detuning,
                    coupling=coupling)
prob_LL = build_problem(comb_LL, DynamicLL(); time_window=1000.0, tpoints=2^12)
sol_LL = solve(prob_LL)
