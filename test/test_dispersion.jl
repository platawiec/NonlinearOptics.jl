#Tests dispersion calculations

f = Frequency.(collect(linspace(c/1000e-9, c/500e-9, 20)))
neff = collect(linspace(1, 2, 20).^2)
aeff = OpticalAttr(1.0e-12, "Aeff")
attr = OpticalAttr(f, neff, "neff")
mode = Mode(attr, aeff)
@test attr(f[10]) ≈ neff[10] # fit function should be close to data
@test NonlinearOptics.der(attr, getω(f[end]); order=0) ≈ 4.0 #derivative at 0 order should
                                             #should return neff
@test get_groupindex(mode, f[1]) ≈ 3.0 # =neff + 2*neff in this special case

attr_nondispersive = OpticalAttr(1.0, "nondispersive")
mode_nondispersive = Mode(attr_nondispersive, aeff)
@test_throws Exception get_beta(mode_nondispersive, f[1], 1)
