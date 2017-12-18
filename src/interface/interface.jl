# This can probably be cleaned up and simplifed with a little thought
@inline wavelength{T<:Unitful.Frequency}(s::T) = c0/s
@inline wavelength{T<:Unitful.Length}(s::T) = s
@inline wavelength(attr::OpticalAttr) = wavelength.(attr.source)
@inline wavelength(laser::AbstractLaser) = wavelength(laser.frequency)
@inline frequency{T<:Unitful.Frequency}(s::T) = s
@inline frequency{T<:Unitful.Length}(s::T) = c0/s
@inline frequency(attr::OpticalAttr) = frequency.(attr.source)
@inline frequency(laser::AbstractLaser) = frequency(laser.frequency)
@inline frequency(model::Model) = frequency(model.laser)
@inline frequency(model::ToyModel) = frequency(model.ω0)

function get_attr(attr::OpticalAttr) attr.property end

function get_label(attr::OpticalAttr) attr.label end
"""
OpticalAttr is callable. Giving a source source will return the
value of the OpticalAttr interpolated at that point
"""
(attr::OpticalAttr)(source) = attr.fit_func(frequency(source))
"""
Alias for der
"""
function der(attr::OpticalAttr, query; order=1)
    der(attr.fit_func, query; order=order)
end

function fit(attr::OpticalAttr, poly_order=12)
    ω = frequency(attr)
    μ = mean(ω)
    σ = std(ω)
    attr_unit = oneunit(first(get_attr(attr)))
    p = attr_unit*polyfit((ω-μ)/σ, get_attr(attr)./attr_unit, poly_order)
    attr.fit_func = ScaledFit(μ, σ, p)
    return attr
end
