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

function get_attr(attr::OpticalAttr) attr.property end

function get_label(attr::OpticalAttr) attr.label end
"""
OpticalAttr is callable. Giving a source source will return the
value of the OpticalAttr interpolated at that point
"""
(attr::OpticalAttr)(source) = attr.fit_func(frequency(source))
(attr::OpticalAttr)(ω::Number) = attr.fit_func(ω)
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

function add_mode!(structure, mode::Mode)
    push!(structure.modes, mode)
end
function add_interaction!(structure, mode_1::Mode, mode_2::Mode, overlap::Float64=1.0)
    !(mode_1 in structure.modes) && error("$mode_1 not found in structure")
    !(mode_2 in structure.modes) && error("$mode_2 not found in structure")

    mode_1_idx = findfirst(structure.modes, mode_1)
    mode_2_idx = findfirst(structure.modes, mode_2)
    mode_1.has_interaction = true
    mode_2.has_interaction = true
    structure.interactions[mode_1_idx] = mode_2_idx
    structure.interactions[mode_2_idx] = mode_1_idx
    structure.overlap[(mode_1_idx, mode_2_idx)] = overlap
    structure.overlap[(mode_2_idx, mode_1_idx)] = overlap
end

"""
    add_straight!(structure, length[, medium])

Adds a straight segment of `length` to `structure`, optionally
changing the medium of the propagating mode to `medium`
"""
function add_straight!(structure, length, medium)
end
function add_straight!(structure, length)
end

function add_turn!(structure, radius, angle)
end

"""
    add_taper!(structure, length)

Adds a tapered segment of `length`. The variables for the parameters
(dispersion, nonlinearity, etc.) in the tapered segment are linearly
interpolated between the parameters for the modes immediately preceding
and succeeding the tapered region. The medium must be the same both
before and after the taper

Not implemented
"""
function add_taper!(structure, length)
    error("add_taper!: Not implemented")
end
