"""
    get_beta -> Vector{Real}

returns the dispersion of the mode at a given wavelength up to
the given order for the structure's modes
"""
function get_beta(effectiveindex, source, numorders::Int)
    ω = getω(effectiveindex)
    length(ω) == 1 && error("the given refractive index is non-dispersive")
    β₀ = ω/c .* get_attr(effectiveindex)
    β = OpticalAttr(frequency(effectiveindex), β₀, "β")
    ω_query = getω(source)
    β_atquery = zeros(eltype(β₀), numorders+1)
    for order=0:numorders
        β_atquery[order+1] = der(β, ω_query; order=order)
    end
    return β_atquery
end

"""
    get_groupindex -> Real

returns the group index of the mode at a given wavelength
"""
function get_groupindex(effindex, source)
    β = get_beta(effindex, source, 1)
    n_group = β[2] * c
    n_group
end

"""
    get_FSR -> Real

returns the free spectral range (FSR) of the resonator's modes
"""
function get_FSR(resonator::AbstractResonator, mode::Mode, source::AbstractSource)
    groupindex = get_groupindex(mode, source)
    FSR = c/(2*groupindex*circumference(resonator))
    FSR
end
