lam = collect(linspace(900e-9, 2000e-9, 250))
neff_data = log10.(collect(linspace(10.0, 50.0, 250)).^2)
neff = OpticalAttr(lam, neff_data, "n_eff")
corefrac = OpticalAttr(1.0, "Core Fraction")
Aeff = OpticalAttr(1e-12, "Aeff (m²)")

mode = Mode(neff, Aeff, corefrac)

model_simple = CircularResonator(20e-6)
add_mode!(model_simple, mode)
