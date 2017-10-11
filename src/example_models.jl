lam = collect(linspace(900e-9, 2000e-9, 250))
neff = log10.(collect(linspace(10.0, 50.0, 250)).^2)
neff_mode = EffectiveRefractiveIndex(lam, neff)
corefrac = collect(linspace(1.0, 1.0, 250))
corefrac_mode = CoreFraction(lam, corefrac)
Aeff = collect(linspace(1.0e-12, 2.0e-12, 250))
Aeff_mode = EffectiveModeArea(lam, Aeff)

mode = Mode(neff_mode, Aeff_mode, corefrac_mode)
