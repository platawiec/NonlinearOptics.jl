#TODO: currently placeholders
const Diamond = Crystal(1.3e-19, 1.0, 1.0)
const Silicon = Crystal(3e-18, 1.0, 1.0)
const SiO2 = Glass(2.7e-20, 1.0, 1.0)

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
    CubicElectronicTensor() -> Tensor{4, 3}

    Returns a dimensionless tensor describing the electronic
    contribution for cubic crystals
"""
function CubicElectronicTensor(ρ)
    δ = (i, j) -> i == j ? 1 : 0
    δ4 = (i, j, k, l) -> i == j && j == k && k == l ? 1 : 0
    f = (i, j, k, l) -> ρ*(δ(i,j)*δ(k,l) + δ(i,k)*δ(j,l)) + (1-ρ)*δ4(i,j,k,l)
    return Tensor{4, 3}(f)
end
