import Base.convert

convert(::Type{Wavelength}, x::Real) = Wavelength(x)
convert(::Type{Wavelength}, f::Frequency) = Wavelength(c/f)
