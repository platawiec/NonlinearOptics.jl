## NLSE Examples

N_linear = (t, u) -> 0.0
D_linear = (k, u) -> 0.0
u0 = z -> (1+0im)*sech.(z)
tspan = (0.0, 1.0)
zmesh = linspace(-100.0, 100.0, 1024)

"""Example of a linear propagating pulse
"""
prob_linear_pulse = NLSEProblem(N_linear, D_linear, u0, tspan, zmesh)


γnl = 9
β₂ = -3

N = (t, u) -> 1im * γnl * abs2.(u) .* u
D = (k, u) -> 1im * β₂ / 2 * k.^2
u0 = z -> (1+0im)*sech.(z)
tspan = (0.0, 1.0)
zmesh = linspace(-100.0, 100.0, 1024)

"""Example of a bright soliton problem construction
"""
prob_bright_soliton = NLSEProblem(N, D, u0, tspan, zmesh)


β₂ = 5e-26
τ = 20e-15
D_dispersive = (k, u) -> 1im * β₂ / 2 * k.^2
u0_gauss = z -> (1+0im)*exp(-z^2/τ^2)
tmesh = linspace(-100e-15, 100e-15, 2^10)
zspan = (0.0, 0.002)
"""A propagating pulse in dispersive fiber example,
with units following SI convention.
"""
prob_dispersive_pulse = NLSEProblem(N_linear, D_dispersive, u0_gauss, zspan, tmesh)

FSR = 100e12
α = 0.01
sqrtθ = sqrt(0.01)
γnl = 100
L = 0.1e-3
Ein = 1
β₂ = -500e-27 #ps²/km -> s²/m
δ₀ = 0.6
D_LL = (k, u) -> 1im*FSR*L/2*k.^2
N_LL = (t, u) -> (1im*γnl*L*FSR*abs2.(u) - (α+1im*δ₀)).*u + sqrtθ*Ein
τspan = (0.0, 1.0)
tmesh = linspace(-500e-15, 500e15, 2^10)
"""Lugatio-Lefever evolution of cavity mean-field. From
"Frequencu Comb Generation beyond the Lugatio-Lefever equation: multi-stability
and super cavity solitons" T. Hansson, S. Wabnitz, Arxiv 1503.03274
Parameters from Fig. 7a
"""
prob_LL = NLSEProblem(N_LL, D_LL, u0_gauss, τspan, tmesh)
