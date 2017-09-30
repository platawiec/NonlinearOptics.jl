#using NonlinearOptics

γnl = 9
β₂ = 3

N = (t, u) -> 1im * γnl * abs2.(u) .* u
D = (k, u) -> 1im * β₂ / 2 * k.^2
u0 = (z) -> sech.(z)
tspan = (0.0, 1.0)
zspan = (-1.0, 1.0)
methods(NLSEProblem)
prob_bright_soliton = NLSEProblem(N, D, u0, tspan, zspan)

# Define a soliton problem
prob = prob_bright_soliton
"""Example of a bright soliton problem construction
"""
# Solve it
println("Split-Step Fourier Method")
sol = solve(prob, SymmetrizedSplitStep())
