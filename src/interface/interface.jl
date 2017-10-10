function wavelength(light::Wavelength) light.λ end
function wavelength(light::Frequency) c/light.f end
function frequency(light::Wavelength) c/light.λ end
function frequency(light::Frequency) light.f end
function frequency(prop::AbstractOpticalProperty) frequency(prop.light) end
function wavelength(prop::AbstractOpticalProperty) wavelength(prop.light) end
function get_property(prop::EffectiveRefractiveIndex) prop.effectiveindex end
function get_property(prop::CoreFraction) prop.corefraction end
function get_property(prop::EffectiveModeArea) prop.effectivearea end

function circumference(res::CircularResonator) 2pi*res.radius end
function circumference(res::RacetrackResonator) 2pi*res.radius + 2*res.length end



"""
    get_beta -> Tuple{Real}

returns the dispersion of the mode at a given wavelength up to
the given order for the structure's modes
"""
function get_beta(mode::Mode, light::AbstractLight, order::Int)

end

"""
    get_groupindex -> Real

returns the group index of the mode at a given wavelength
"""
function get_groupindex(mode::Mode, light::AbstractLight)
end

"""
    get_FSR -> Real

returns the free spectral range (FSR) of the resonator's modes
"""
function get_FSR(resonator::AbstractResonator, mode::Mode, light::AbstractLight)
    groupindex = get_groupindex(mode, light)
    FSR = c/(2*groupindex*circumference(resonator))
    FSR
end
