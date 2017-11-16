# This can probably be cleaned up and simplifed with a little thought
function wavelength(source::Wavelength) source.λ end
function wavelength(source::Frequency) c/source.f end
function wavelength(attr::AbstractOpticalAttr) wavelength.(attr.source) end
function wavelength(x) wavelength(convert(Wavelength, x)) end
function wavelength(laser::AbstractLaser) wavelength(laser.frequency) end
function frequency(source::Wavelength) c/source.λ end
function frequency(source::Frequency) source.f end
function frequency(attr::AbstractOpticalAttr) frequency.(attr.source) end
function frequency(laser::AbstractLaser) frequency(laser.frequency) end
function frequency(x) frequency(convert(Wavelength, x)) end
function getω(source) 2pi*frequency(source) end

function get_attr(attr::OpticalAttr) attr.property end

function get_label(attr::OpticalAttr) attr.label end
"""
OpticalAttr is callable. Giving a source source will return the
value of the AbstractOpticalAttr interpolated at that point
"""
(attr::OpticalAttr)(source) = attr.fit_func(getω(source))
"""
Alias for der
"""
function der(attr::AbstractOpticalAttr, query; order=1)
    der(attr.fit_func, query; order=order)
end

function fit(attr::AbstractOpticalAttr, poly_order=12)
    ω = getω(attr)
    μ = mean(ω)
    σ = std(ω)
    p = polyfit((ω-μ)/σ, get_attr(attr), poly_order)
    attr.fit_func = ScaledFit(μ, σ, p)
    attr
end

function add_mode!(structure, mode::Mode)
    push!(structure.modes, mode)
end
