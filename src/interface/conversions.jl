import Base.convert

convert(::Type{AbstractSource}, x::Real) = Frequency(x)
convert(::Type{Wavelength}, x::Real) = Wavelength(x)
convert(::Type{Wavelength}, f::Frequency) = Wavelength(c/f.f)
convert(::Type{Frequency}, x::Real) = Frequency(x)
convert(::Type{Frequency}, λ::Wavelength) = Frequency(c/λ.λ)
