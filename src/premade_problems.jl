## NLSE Examples
γnl = 9
β₂ = 3

function N(t, u) 1im * γnl * abs2.(u) .* u end
function D(k, u) 1im * β₂ / 2 * k.^2 end
u0 = sech.(-1.0:0.01:1.0)
tspan = (0.0, 1.0)
zspan = (-1.0, 1.0)

"""Example of a bright soliton problem construction
"""
prob_bright_soliton = NLSEProblem(N, D, u0, tspan, zspan)
