# Tests for the ways various types interface with eachother

@test frequency(1.0) ≈ c
@test wavelength(1.0) ≈ 1.0
@test wavelength(BigFloat(1.0)) ≈ 1.0

attr = OpticalAttr(1.0, "Attribute")
@test attr(1.0) == attr(500e-9)
@test attr(Wavelength(200e-9)) == attr(Frequency(150e12))
