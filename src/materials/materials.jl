"""
    CubicRamanTensor() -> Tensor{4, 3}

    Returns a dimensionless tensor describing the Raman
    contribution for cubic crystals
"""
function CubicRamanTensor()
    δ = (i, j) -> i == j ? 1 : 0
    δ4 = (i, j, k, l) -> i == j && j == k && k == l ? 1 : 0
    f = (i, j, k, l) -> δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k) - 2*δ4(i,j,k,l)
    return Tensor{4, 3}(f)
end

"""
    CubicElectronicTensor(ρ) -> Tensor{4, 3}

    Returns a dimensionless tensor describing the electronic
    contribution for cubic crystals. ρ characterizes the nonlinear
    anisotropy ρ ≡ 3χ₁₁₂₂/χ₁₁₁₁
"""
function CubicElectronicTensor(ρ)
    δ = (i, j) -> i == j ? 1 : 0
    δ4 = (i, j, k, l) -> i == j && j == k && k == l ? 1 : 0
    f = (i, j, k, l) -> ρ/3*(δ(i,j)*δ(k,l) + δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + (1-ρ)*δ4(i,j,k,l)
    return Tensor{4, 3}(f)
end

"""
    IsotropicTensor(λ, μ, ν) -> Tensor{4, 3}
    IsotropicTensor()        -> Tensor{4, 3}

    Returns a dimensionless tensor describing the
    an isotropic nonlinear contribution. λ, μ, and ν are
    parameters which control the coefficients, see
    https://farside.ph.utexas.edu/teaching/336L/Fluid/node252.html
"""
function IsotropicTensor(λ, μ, ν)
    δ = (i, j) -> i == j ? 1 : 0
    f = (i, j, k, l) -> (λ*δ(i,j)*δ(k,l) + μ*δ(i,k)*δ(j,l) + ν*δ(i,l)*δ(j,k))
    return Tensor{4, 3}(f)
end
IsotropicTensor() = IsotropicTensor(1,1,1)

const Diamond = Material(
                    ElectronicTensor(
                        CubicElectronicTensor(1.2),
                        1.3e-19m^2/W
                    ),
                    RamanTensor(
                        CubicRamanTensor(),
                        0.004ps,
                        5.7ps,
                        0.28
                    )
                )
const Silicon = Material(
                    ElectronicTensor(
                        CubicElectronicTensor(1.27),
                        3.8e-18m^2/W
                    ),
                    RamanTensor(
                        CubicRamanTensor(),
                        0.01ps,
                        3.0ps,
                        0.08
                    )
                )
const SiO2 = Material(
                    ElectronicTensor(
                        IsotropicTensor(),
                        2.7e-20m^2/W,
                    ),
                    RamanTensor(
                        IsotropicTensor(),
                        0.0122ps,
                        0.032ps,
                        0.18
                    )
                )
