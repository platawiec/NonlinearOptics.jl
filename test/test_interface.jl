# Tests for the ways various types interface with eachother

# Test conversions
@test frequency(1.0m)*m ≈ c0
@test frequency(1nm)*nm ≈ c0
@test frequency(wavelength(100THz)) ≈ 100THz
@test wavelength(frequency(100nm)) ≈ 100nm
@test wavelength(100nm) ≈ wavelength(100e-9m)

attr = OpticalAttr(1.0, "Attribute")
@test attr(1.0THz) == attr(500nm)
@test attr(100THz) == attr(200THz)
@test attr(frequency(200nm)) == attr(frequency(150THz))
@test NonlinearOptics.der(attr, 100nm) == 0.0
@test NonlinearOptics.der(attr, 100nm; order=0) == 1.0
@test NonlinearOptics.der(attr, 1.0nm; order=3) == 0.0
@test NonlinearOptics.der(attr, 100THz; order=3) == 0.0

frequencies = linspace(100THz, 200THz, 20)
unitful_attr = frequencies./THz.*m
unitful_attr = OpticalAttr(frequencies, unitful_attr, "Attribute [m]")
@test unitful_attr(150THz) ≈ 150m
@test unitful_attr(frequency(c0/150THz)) ≈ 150m
@test unitful_attr(frequency(1998.6163866nm)) ≈ 150m
@test unitful_attr(1998.6163866nm) ≈ 150m
@test NonlinearOptics.der(unitful_attr, 150THz, order=0) == unitful_attr(150THz)
@test NonlinearOptics.der(unitful_attr, 150THz, order=1) ≈ 1m/THz
@test isapprox(NonlinearOptics.der(unitful_attr, 150THz, order=2), 0.0m/THz^2; atol=5e-15m/THz^2)
