# Tests for the ways various types interface with eachother

@test frequency(1.0) ≈ c
@test wavelength(1.0) ≈ 1.0
@test wavelength(BigFloat(1.0)) ≈ 1.0
