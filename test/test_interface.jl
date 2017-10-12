# Tests for the ways various types interface with eachother

@test frequency(1.0) ≈ c
@test wavelength(1.0) ≈ 1.0
@test wavelength(BigFloat(1.0)) ≈ 1.0

attr = OpticalAttr(1.0, "Attribute")
@test attr(1.0) == attr(500e-9)
@test attr(Wavelength(200e-9)) == attr(Frequency(150e12))
@test NonlinearOptics.der(attr, 100e-9) == 0.0
@test NonlinearOptics.der(attr, 100e-9; order=0) == 1.0
@test NonlinearOptics.der(attr, 1.0; order=3) == 0.0

f = Frequency.(collect(linspace(c/1000e-9, c/500e-9, 20)))
neff = collect(linspace(1, 2, 20).^2)
attr = OpticalAttr(f, neff, "neff")
@test attr(f[10]) ≈ neff[10] # fit function should be close to data
@test der(attr, getω(f[end]); order=0) ≈ 4.0 #derivative at 0 order should
                                             #should return neff
@test get_groupindex(attr, f[1]) ≈ 3.0 # =neff + 2*neff in this case

get_beta(attr, f[10], 4)
attr(f[10])
der(attr, getω(f[end]); order=1)
attr(f[10])
get_beta(attr, f[end], 0)
get_groupindex(attr, f[end])
