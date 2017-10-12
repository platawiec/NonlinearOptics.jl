# Tests for the ways various types interface with eachother

@test frequency(1.0) ≈ c
@test wavelength(1.0) ≈ 1.0
@test wavelength(BigFloat(1.0)) ≈ 1.0

attr = OpticalAttr(1.0, "Attribute")
@test attr(1.0) == attr(500e-9)
@test attr(Wavelength(200e-9)) == attr(Frequency(150e12))
@test NonlinearOptics.der(attr, 100e-9) == 0
@test NonlinearOptics.der(attr, 100e-9; order=0) == 1.0
@test NonlinearOptics.der(attr, 1.0; order=3) == 0
