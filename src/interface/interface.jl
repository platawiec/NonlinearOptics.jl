function wavelength(source::Wavelength) source.λ end
function wavelength(source::Frequency) c/source.f end
function frequency(source::Wavelength) c/source.λ end
function frequency(source::Frequency) source.f end
function frequency(prop::AbstractOpticalProperty) frequency.(prop.source) end
function wavelength(prop::AbstractOpticalProperty) wavelength.(prop.source) end
function get_property(prop::EffectiveRefractiveIndex) prop.effectiveindex end
function get_property(prop::CoreFraction) prop.corefraction end
function get_property(prop::EffectiveModeArea) prop.effectivearea end
function get_label(prop::EffectiveModeArea) return "Effective Mode Area (m²)" end
function get_label(prop::CoreFraction) return "Core Fraction" end
function get_label(prop::EffectiveRefractiveIndex) "Effective Refractive Index" end
"""
AbstractOpticalProperty is callable. Giving a source source will return the
value of the AbstractOpticalProperty interpolated at that point
"""
function (prop::AbstractOpticalProperty)(source::AbstractSource)
    prop.poly_interp(2pi*frequency(source))
end


function circumference(res::CircularResonator) 2pi*res.radius end
function circumference(res::RacetrackResonator) 2pi*res.radius + 2*res.length end



"""
    get_beta -> Vector{Real}

returns the dispersion of the mode at a given wavelength up to
the given order for the structure's modes
"""
function get_beta(mode::Mode, source::AbstractSource, numorders::Int)
    ω = 2pi * frequency(mode)
    β₀ = ω/c .* get_property(mode.effectiveindex)
    σ = std(ω)
    μ = mean(ω)
    ω_rescaled = (ω-μ)/σ
    p = polyfit(ω_rescaled, β₀)
    ω_query = (2pi*(frequency(source)) - μ)/σ
    β_atquery = zeros(typeof(β₀), numorders)
    for i=0:numorders
        β_atquery[i+1] = polyder(p, order=i)(ω_query)/(σ^order)
    end
    return β_atquery
end

"""
    get_groupindex -> Real

returns the group index of the mode at a given wavelength
"""
function get_groupindex(mode::Mode, source::AbstractSource)
    β = get_beta(mode, source, 1)
    n_group = β[1] + 2pi*frequency(source)*β[2]
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
