## NLSE Examples
f(t, z) -> sech.(z)
A₀ = z -> sech.(z)
analytic_soliton(t, z) = sech.(z)
tspan = (0.0, 1.0)
zspan = (0.0, 1.0)
dz = 1e-5
dt = 1e-5

"""Example problem for a bright soliton (stationary solution):

```math
A(t, z)=sech(z)
```
"""
prob_bright_soliton = NLSEProblem(analytic_soliton, A₀)
