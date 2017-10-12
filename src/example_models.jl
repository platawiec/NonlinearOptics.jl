lam = collect(linspace(900e-9, 2000e-9, 250))
neff = log10.(collect(linspace(10.0, 50.0, 250)).^2)
neff_mode = OpticalAttr(lam, neff, "n_eff")
corefrac = collect(linspace(1.0, 1.0, 250))
corefrac_mode = OpticalAttr(lam, corefrac, "Core Fraction")
Aeff = collect(linspace(1.0e-12, 2.0e-12, 250))
Aeff_mode = OpticalAttr(lam, Aeff, "Aeff (m²)")

mode = Mode(neff_mode, Aeff_mode, corefrac_mode)

#get_groupindex(mode, Wavelength(700e-9))
