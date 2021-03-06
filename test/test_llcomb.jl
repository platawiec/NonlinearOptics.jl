# Simple LL comb simulation
T_FWHM = 1.0ps
T_T0 = T_FWHM/(2*log(1+sqrt(2)))
β₂ = -0.013ps^2/m
γnl = 0.000032W/m
α = 1.75e-5/m
coupling = 1.75e-5
P_in= 55.6e-3W
FSR = 0.0182THz
L = 0.0119m
detuning = 0.0012
comb_LL = ToyModel(;betacoeff=[0.0/m, 0.0ps/m, β₂],
                    nonlinearcoeff=γnl,
                    linearloss=α,
                    power_in=P_in,
                    pulsetime=0.0ps,
                    length=L,
                    FSR=FSR,
                    detuning=detuning,
                    coupling=coupling)
prob_LL = build_problem(comb_LL, DynamicLL(); time_window=10.0ps, tpoints=2^12)
sol_LL = solve(prob_LL)
