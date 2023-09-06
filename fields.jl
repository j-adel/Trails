include("fns.jl")
include("fieldFns.jl")

function vectorfield(p)
    # Initialize the arrays to store s and w values
    vectors = Point[]
    weights = []

    # s, w = attractorLine(p, Point(0, -W/4), Point(0, W/4), 0.7, 150)
    # push!(vectors, s)
    # push!(weights, w)

    # s, w = attractorLine(p, Point(-W/4, 0), Point(W/4, 0), -0.7, 150)
    # push!(vectors, s)
    # push!(weights, w)

    s, w = attractorPoint(p, Point(0, 0), 0.9, 150)
    push!(vectors, s)
    push!(weights, w)

    # Compute the weighted average of the s values
    weighted_sum_s = sum(weights[i] * vectors[i] for i in eachindex(vectors))
    total_w = sum(weights) + 0.01
    s_avg = weighted_sum_s / total_w

    s_avg /= (W/2)
    s_avg *= 15
    return s_avg
end


abstract type Field end

mutable struct attPointField <: Field
    f::Function
    c::Point
    rotation::Float64
    areaFactor::Float64
    seedPoints::Vector{Point}
end

mutable struct attLineField <: Field
    f::Function
    A::Point
    B::Point
    rotation::Float64
    areaFactor::Float64
    seedPoints::Vector{Point}
end

    
