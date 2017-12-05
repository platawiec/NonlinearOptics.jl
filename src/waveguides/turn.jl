"""
    struct Turn{T} <: Segment{T}
        axis::Tensors.Vec{3, Float64}
        α::typeof(1.0°)
        r::T
        p0::T
        α0::typeof(1.0°)
    end
A circular turn is parameterized by the turn angle `α` and turning radius `r`.
It begins initial angle `α0` around the axis `axis`.
The center of the circle is given by:
`cen = p0 + Point(r*cos(α0+sign(α)*π/2), r*sin(α0+sign(α)*π/2))`
The parametric function over `t ∈ [0,1]` describing the turn is given by:
`t -> cen + Point(r*cos(α0-sign(α)*π/2+α*t), r*sin(α0-sign(α)*π/2+α*t))`
"""
struct Turn{T} <: Segment{T}
    α::typeof(1.0°)
    r::T
    p0::T
    axis::Tensors.Vec{3, Float64}
    α0::typeof(1.0°)
end

"""
    Turn{T<:Unitful.Length}(α, r::T, p0::T,
        axis=Tensors.Vec{3}((0., 0., 1.)), α0=0.0°)
Outer constructor for `Turn` segments.
"""
Turn(α, r::T;
     p0::T=zero(T),
     axis=Tensors.Vec{3}((0., 0., 1.)),
     α0=0.0°) where {T <: Unitful.Length} =
    Turn{T}(α, r, p0, axis, α0)

pathlength(s::Turn{T}) where {T} = T(abs(s.r*s.α))
function orientation(s::Turn, z)
    axis = s.axis
    angle = s.α0 + (z-s.p0)/s.r
    return (axis, angle)
end

"""
    add_turn!(w::Waveguide{T}, α, r::Unitful.Length,
        med::Medium=lastmedium(w)) where {T <: Unitful.Length}
Turn a waveguide `w` by angle `α` with turning radius `r` in the
current direction. By default, we take the last continuous medium
in the waveguide
"""
function add_turn!(w::Waveguide{T}, α, r::Unitful.Length,
        med::Medium=lastmedium(w)) where {T <: Unitful.Length}
    Unitful.dimension(T) != Unitful.dimension(typeof(r)) && throw(Unitful.DimensionError(T(1),l))
    @assert r > zero(r) "tried to turn by infinitesimal amount."
    p0 = pathlength(w)
    axis, α0 = orientation(w, p0)
    s = Turn{T}(α, r, p0, axis, α0)
    push!(w.structures, Structure(s, med))
    return nothing
end
