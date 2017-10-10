function wavelength(light::Wavelength) light.λ end
function wavelength(light::Frequency) c/light.f end
function frequency(light::Wavelength) c/light.λ end
function frequency(light::Frequency) light.f end

"""
    get_beta -> Tuple{Real}

returns the dispersion of the structure at a given wavelength up to
the given order
"""
function get_beta(structure::AbstractStructure, light::AbstractLight, order::Int)

end
