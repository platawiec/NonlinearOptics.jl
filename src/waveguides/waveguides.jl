abstract type Segment{T<:Unitful.Length} end
@inline Base.eltype(::Segment{T}) where {T} = T
@inline Base.eltype(::Type{Segment{T}}) where {T} = T

struct NullSegment{T} <: Segment{T}
    p0::T
    axis::Tensors.Vec{3, Float64}
    α0::typeof(1.0°)
end
NullSegment(p0::T;
            axis=Tensors.Vec{3}((0.0, 0.0, 1.0)),
            α0=0.0°) where {T<:Unitful.Length} = NullSegment{T}(p0, axis, α0)
NullSegment() = NullSegment(0.0m)

pathlength(::NullSegment{T}) where T = zero(T)
orientation(ns::NullSegment, z) = (ns.axis, ns.α0)

struct Medium
    material::Material
    modes::Vector{Mode}
    interactions::Dict{Int, Int}
    overlap::Dict{Tuple{Int, Int}, Float64}
end
Medium(material) = Medium(material, Mode[],
                       Dict{Int, Int}(),
                       Dict{Tuple{Int, Int}, Float64}()
                   )

struct Structure{T<:Unitful.Length}
    segment::Segment{T}
    medium::Medium
end

struct Waveguide{T<:Unitful.Length}
    structures::Vector{Structure}
end

"""
    Waveguide()

Outer constructors for `Waveguide`
"""
Waveguide(material::Material) = Waveguide{eltype(0.0m)}(
                                    [Structure(
                                        NullSegment(),
                                        Medium(material)
                                    )]
                                )

pathlength(w::Waveguide) = sum(pathlength, w.structures)
pathlength(s::Structure) = pathlength(s.segment)
function currentstructure(w::Waveguide, z)
    z_init = zero(z)
    idx = 1
    last_struct = last(w.structures)
    for s in w.structures
        if z_init + pathlength(s) >= z
            return idx, s
        end
        idx    += 1
        z_init += pathlength(s)
    end
    return idx, last_struct
end
function orientation(w::Waveguide, z)
    idx, structure = currentstructure(w, z)
    return orientation(structure, z)
end
orientation(s::Structure, z) = orientation(s.segment, z)

lastmedium(w::Waveguide) = last(w.structures).medium
currentmedium(w::Waveguide, z) = currentstructure(w, z)[2].medium
newmedium(m::Material) = Medium(m)

function add_mode!(w::Waveguide, mode::Mode)
    add_mode!(lastmedium(w), mode)
end
add_mode!(m::Medium, mode::Mode) = push!(m.modes, mode)

function add_interaction!(w::Waveguide, mode_1::Mode, mode_2::Mode, overlap::Float64=1.0)
    add_interaction!(lastmedium(w), mode_1, mode_2, overlap)
end
function add_interaction!(med::Medium, mode_1::Mode, mode_2::Mode, overlap::Float64=1.0)
    !(mode_1 in med.modes) && error("$mode_1 not found in medium")
    !(mode_2 in med.modes) && error("$mode_2 not found in medium")

    mode_1_idx = findfirst(med.modes, mode_1)
    mode_2_idx = findfirst(med.modes, mode_2)
    mode_1.has_interaction = true
    mode_2.has_interaction = true
    med.interactions[mode_1_idx] = mode_2_idx
    med.interactions[mode_2_idx] = mode_1_idx
    med.overlap[(mode_1_idx, mode_2_idx)] = overlap
    med.overlap[(mode_2_idx, mode_1_idx)] = overlap
end

@inline num_modes(w::Waveguide) = num_modes(first(w.structures))
@inline num_modes(s::Structure) = num_modes(s.medium)
@inline num_modes(m::Medium) = length(m.modes)
@inline num_structures(m::Waveguide) = length(m.structures)

include("straight.jl")
include("turn.jl")
include("taper.jl")
