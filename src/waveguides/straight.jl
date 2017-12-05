"""
    struct Straight{T} <: ContinuousSegment{T}
        l::T
        α0::typeof(0.0°)
        axis::Tensors.Vec{3, Float64}
    end
A straight line segment is parameterized by its length.
It begins at a point `p0` with initial angle `α0`.
The parametric function describing the line segment is given by
`t -> p0 + Point(t*cos(α),t*sin(α))` where `t` is a length from 0 to `l`.
"""
struct Straight{T} <: Segment{T}
    l::T
    p0::T
    axis::Tensors.Vec{3, Float64}
    α0::typeof(1.0°)
end

"""
    Straight{T<:Coordinate}(l::T; p0::Point=Point(zero(T),zero(T)), α0=0.0°)
Outer constructor for `Straight` segments.
"""
Straight(l::T;
         p0::T=zero(T),
         axis=Tensors.Vec{3}((0.0, 0.0, 1.0)),
         α0=0.0°) where {T <: Unitful.Length} =
    Straight{T}(l, p0, axis, α0)

pathlength(s::Straight) = s.l
orientation(s::Straight, z) = (s.axis, s.α0)

"""
    add_straight!(w::Waveguide{T}, l::Unitful.Length,
        med::Medium=lastmedium(w)) where {T <: Unitful.Length}
Extend a waveguide `w` straight by length `l` in the current direction.
By default, we take the last continuous medium in the waveguide
"""
function add_straight!(w::Waveguide{T}, l::Unitful.Length,
        med::Medium=lastmedium(w)) where {T <: Unitful.Length}
    Unitful.dimension(T) != Unitful.dimension(typeof(l)) && throw(Unitful.DimensionError(T(1),l))
    @assert l > zero(l) "tried to go straight by zero or a negative amount."
    p0 = pathlength(w)
    (axis, α) = orientation(w, p0)
    s = Straight{T}(l, p0, axis, α)
    push!(w.structures, Structure(s, med))
    return nothing
end
