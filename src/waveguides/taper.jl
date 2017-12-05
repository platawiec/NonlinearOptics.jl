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
