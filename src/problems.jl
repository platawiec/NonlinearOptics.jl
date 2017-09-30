"""
∂A/∂t = N A + D A
D is dispersive operator, N is nonlinear operator
Note that for NLSE for pulse propagation, substitute z -> t
"""

mutable struct NLSEProblem{uType, tType, kType, F1, F2, C} <: AbstractNLSEProblem{uType, tType, zType}
    N::F1
    D::F2
    u0::uType
    tspan::Tuple{tType,tType}
    zspan::Tuple{zType,zType}
    callback::C
    function NLSEProblem(N, D, u0, tspan, zspan; callback=nothing)
        new{typeof(u0), promote_type(map(typeof, tspan)...),
            promote_type(map(typeof, zspan)...),
            typeof(N), typeof(D)}(
            N, D, u0, tspan, zspan, callback)
    end
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
