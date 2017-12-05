"""
    get_beta -> Vector{Real}

returns the dispersion of the mode at a given wavelength up to
the given order for the structure's modes
"""
function get_beta(mode::Mode, source)
    effectiveindex = mode.effectiveindex
    ω = frequency(effectiveindex)
    length(ω) == 1 && error("the given refractive index is non-dispersive")
    β₀ = ω/c0 .* get_attr(effectiveindex)
    β = OpticalAttr(ω, β₀, "β")
    return β
end

function get_beta(mode::Mode, source, numorders::Int)
    effectiveindex = mode.effectiveindex
    ω = frequency(effectiveindex)
    length(ω) == 1 && error("the given refractive index is non-dispersive")
    β₀ = ω/c0 .* get_attr(effectiveindex)
    β = OpticalAttr(frequency(effectiveindex), β₀, "β")
    ω_query = frequency(source)
    β_atquery = Any[]
    for order=0:numorders
        push!(β_atquery, der(β, ω_query; order=order))
    end
    return β_atquery
end

function get_beta(mode::ToyMode, source, numorders::Int)
    β = mode.beta
    (numorders+1) > length(β) && error("model does not have sufficient orders of β")
    return β[1:(numorders+1)]
end

"""
    get_stationarydispersion(mode, source)

    Returns the stationary dispersion of the mode at the frequency defined by
    source. This dispersion is centered around ω=0 and has no constant or
    linear component - it is the dispersion in the moving frame.
"""
function get_stationarydispersion(mode, source)
    const source_ω = frequency(source)
    const beta = get_beta(mode, source)
    const beta_linear = get_beta(mode, source, 1)
    return ω->(beta(ω) - beta_linear[1] - beta_linear[2] * (ω - source_ω))
end

function get_dispersionfrombeta(beta_coeff)
    return Poly(beta_coeff, :ω)
end

"""
    get_groupindex -> Real

returns the group index of the mode at a given wavelength
"""
function get_groupindex(mode, source)
    β = get_beta(mode, source, 1)
    n_group = β[2] * c0
    n_group
end

"""
    get_FSR -> Real

returns the free spectral range (FSR) of the resonator's modes
"""
function get_FSR(resonator, mode::Mode, source)
    groupindex = get_groupindex(mode, source)
    FSR = abs(c0/(groupindex*pathlength(resonator)))
    FSR
end
