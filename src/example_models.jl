lam = collect(linspace(900e-9, 2000e-9, 250))
neff_data = log10.(collect(linspace(10.0, 50.0, 250)).^2)
neff = OpticalAttr(lam, neff_data, "n_eff")
corefrac = OpticalAttr(1.0, "Core Fraction")
Aeff = OpticalAttr(1e-12, "Aeff (m²)")
loss = OpticalAttr(1e-3, "Loss (1/m)")

mode = Mode(neff, Aeff, loss, corefrac)

res_simple = CircularResonator(20e-6)
add_mode!(res_simple, mode)
#source_simple = Laser()
#model_simple = Model(res_simple, source_simple)
# Laser answers: δ₀ can change during experiment, ω₀, and Ein? All can be
# parameters of cw laser. Solver ?s: How many dispersive orders
# should we use (keyword, default=3)? Initial conditions for NLSE? How do we
# determine initial conditions otherwise?
#solve(model_simple, DynamicLL()) (DynamicLL, DynamicIkeda, SSLL, SSIkeda, NLSE)
