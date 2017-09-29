"""
∂A/∂z = -𝚤β₂/2 ∂²A/∂t² + 𝚤γ|A|²A
"""
mutable struct NLSEProblem{uType, tType, zType} <: AbstractNLSEProblem{uType, tType, zType}
    u0::uType
    tspan::Tuple{tType,tType}
    zspan::Tuple{zType,zType}
end

"""
Full dispersion
∂A/∂z = ∑₂ⁿ(𝚤ⁱ𝚤 βᵢ/i! ∂ⁱA/∂tⁱ) - α/2 A + 𝚤γ|A|²A
"""

"""
Full dispersion and basic nonlinear terms
∂A/∂z = ∑₂ⁿ(𝚤ⁱ𝚤 βᵢ/i! ∂ⁱA/∂tⁱ) - α/2 A
        + 𝚤γ[|A|²A + 𝚤/ω₀∂(|A|²A)/∂t - Tᵣ∂(|A|²)/∂t A]
"""

"""
Full dispersion and Raman response
∂A/∂z = ∑₂ⁿ(𝚤ⁱ𝚤 βᵢ/i! ∂ⁱA/∂tⁱ) - α/2 A 
        + 𝚤γ[|A|²A + 𝚤/ω₀∂(|A|²A)/∂t - Tᵣ∂(|A|²)/∂t A]
"""
