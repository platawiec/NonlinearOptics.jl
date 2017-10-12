function wavelength(source::Wavelength) source.λ end
function wavelength(source::Frequency) c/source.f end
function frequency(source::Wavelength) c/source.λ end
function frequency(source::Frequency) source.f end
function frequency(attr::AbstractOpticalAttr) frequency.(attr.source) end
function wavelength(attr::AbstractOpticalAttr) wavelength.(attr.source) end
function wavelength(x) wavelength(convert(Wavelength, x)) end
function frequency(x) frequency(convert(Wavelength, x)) end
function getω(source) 2pi*frequency(source) end

function get_attr(attr::OpticalAttr) attr.property end

function get_label(attr::OpticalAttr) attr.label end
"""
OpticalAttr is callable. Giving a source source will return the
value of the AbstractOpticalAttr interpolated at that point
"""
(attr::OpticalAttr)(source::AbstractSource) = attr.fit_func(getω(source))
"""
Alias for der
"""
function der(attr::AbstractOpticalAttr, query; order=1)
    der(attr.fit_func, query; order=order)
end


function circumference(res::CircularResonator) 2pi*res.radius end
function circumference(res::RacetrackResonator) 2pi*res.radius + 2*res.length end

function fit(attr::AbstractOpticalAttr, poly_order=12)
    ω = getω(attr)
    μ = mean(ω)
    σ = std(ω)
    p = polyfit((ω-μ)/σ, get_attr(attr), poly_order)
    attr.fit_func = ScaledFit(μ, σ, p)
    attr
end

"""
    get_beta -> Vector{Real}

returns the dispersion of the mode at a given wavelength up to
the given order for the structure's modes
"""
function get_beta(mode::Mode, source::AbstractSource, numorders::Int)
    ω = getω(mode.effectiveindex)
    β₀ = ω/c .* get_attr(mode.effectiveindex)
    β = OpticalAttr(frequency(mode.effectiveindex), β₀, "β")
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
function get_groupindex(mode::Mode, source::AbstractSource)
    β = get_beta(mode, source, 1)
    n_group = β[1] + getω(source)*β[2]
    return n_group
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
