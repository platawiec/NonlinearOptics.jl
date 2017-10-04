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
tmesh = linspace(-100e-15, 100e-15, 2^12)
zspan = (0.0, 0.001)
prob_dispersive_pulse = NLSEProblem(N_linear, D_dispersive, u0_gauss, zspan, tmesh)
