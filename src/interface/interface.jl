function wavelength(source::Wavelength) source.λ end
function wavelength(source::Frequency) c/source.f end
function frequency(source::Wavelength) c/source.λ end
function frequency(source::Frequency) source.f end
function frequency(prop::AbstractOpticalProperty{T}) frequency.(prop.source) end
function wavelength(prop::AbstractOpticalProperty{T}) wavelength.(prop.source) end
function getω(source) 2pi*frequency(source) end

function get_property(prop::EffectiveRefractiveIndex) prop.effectiveindex end
function get_property(prop::CoreFraction) prop.corefraction end
function get_property(prop::EffectiveModeArea) prop.effectivearea end
function get_property(prop::GenericOpticalProperty) prop.property end

function get_label(prop::EffectiveModeArea) "Effective Mode Area (m²)" end
function get_label(prop::CoreFraction) "Core Fraction" end
function get_label(prop::EffectiveRefractiveIndex) "Effective Refractive Index" end
function get_label(prop::GenericOpticalProperty) prop.label end
"""
AbstractOpticalProperty is callable. Giving a source source will return the
value of the AbstractOpticalProperty interpolated at that point
"""
function (prop::AbstractOpticalProperty{T})(source::AbstractSource)
    prop.fit_func(getω(source))
end
"""
Alias for der
"""
function der(prop::AbstractOpticalProperty{T}, query, order=1)
    der(prop.fit_func, query, order=order)
end


function circumference(res::CircularResonator) 2pi*res.radius end
function circumference(res::RacetrackResonator) 2pi*res.radius + 2*res.length end

function fit(prop::AbstractOpticalProperty{T})
    ω = getω(s)
    μ = mean(ω)
    σ = std(ω)
    p = polyfit((ω-μ)/σ, get_property(prop))
    ScaledFit(μ, σ, p)
end

"""
    get_beta -> Vector{Real}

returns the dispersion of the mode at a given wavelength up to
the given order for the structure's modes
"""
function get_beta(mode::Mode, source::AbstractSource, numorders::Int)
    β₀ = ω/c .* get_property(mode.effectiveindex)
    β = GenericOpticalProperty(mode.f, β₀, "β")
    ω_query = getω(source)
    β_atquery = zeros(typeof(β₀), numorders)
    for order=0:numorders
        β_atquery[order+1] = der(β, order=order)(ω_query)
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
