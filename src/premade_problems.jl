## NLSE Examples
γ = 9
β₂ = 3

N = (t, u) -> 1im * γ * abs2.(u) .* u
D = (k, u) -> 1im * β₂ / 2 * k.^2
u0 = (z) -> sech.(z)
tspan = (0, 1.0)
zspan = (-1.0, 1.0)

"""Example problem for a bright soliton (stationary solution):

```math
A(t, z)=sech(z)
```
"""
prob_bright_soliton = NLSEProblem(N, D, tspan, zspan)
