#TODO: currently placeholders
const Diamond = Crystal(1.3e-19, 0.004, 5.7)
const Silicon = Crystal(3e-18, 0.01, 3.0)
const SiO2 = Glass(2.7e-20, 0.0122, 0.032)

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
    f = (i, j, k, l) -> ρ*(δ(i,j)*δ(k,l) + δ(i,k)*δ(j,l)) + (1-ρ)*δ4(i,j,k,l)
    return Tensor{4, 3}(f)
end
