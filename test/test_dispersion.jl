#Tests dispersion calculations

f = linspace(100THz, 200THz, 20)
neff = collect(linspace(1, 2, 20).^2)
aeff = OpticalAttr(1.0e-12m^-2, "Aeff")
attr = OpticalAttr(f, neff, "neff")
mode = Mode(attr, aeff)
@test attr(f[10]) ≈ neff[10] # fit function should be close to data
@test NonlinearOptics.der(attr, f[end]; order=0) ≈ 4.0
# =neff + 2*neff in this special case
@test get_groupindex(mode, f[1]) ≈ 3.0
attr_nondispersive = OpticalAttr(1.0, "nondispersive")
mode_nondispersive = Mode(attr_nondispersive, aeff)
@test_throws Exception get_beta(mode_nondispersive, f[1], 1)
