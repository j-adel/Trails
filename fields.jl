include("fns.jl")
abstract type Field end

mutable struct attPointField <: Field
    fieldFunction::Function
    c::Point
    rotation::Float64
    areaFactor::Float64
    seedPoints::Vector{Point}
end

mutable struct attLineField <: Field
    fieldFunction::Function
    A::Point
    B::Point
    rotation::Float64
    areaFactor::Float64
    seedPoints::Vector{Point}
end

# function vectorfield(p)
#     # Initialize the arrays to store s and w values
#     vectors = Point[]
#     weights = []

#     s, w = attractorLine(p, Point(0, -W/4), Point(0, W/4), 0.7, 50)
#     push!(vectors, s)
#     push!(weights, w)

#     s, w = attractorLine(p, Point(-W/4, 0), Point(W/4, 0), -0.7, 50)
#     push!(vectors, s)
#     push!(weights, w)

#     # s, w = attractorPoint(p, Point(0, 0), 0.5, 150)
#     # push!(vectors, s)
#     # push!(weights, w)

#     # Compute the weighted average of the s values
#     weighted_sum_s = sum(weights[i] * vectors[i] for i in eachindex(vectors))
#     total_w = sum(weights) + 0.01
#     s_avg = weighted_sum_s / total_w

#     s_avg /= (W/2)
#     s_avg *= 20
#     return s_avg
# end

function initializeFields!(fields)
    field1=attLineField(attractorLine, Point(0, -W/4), Point(0, W/4), 0.7, 50, Point[])
    push!(fields, field1)
    field2=attLineField(attractorLine, Point(-W/4, 0), Point(W/4, 0), -0.7, 50, Point[])
    push!(fields, field2)
    # push!(fields, attPointField((p)->vectorfield(p), Point(0, 0), 0.5, 150, Point[]))
end

function computeVector(p,fields)
    vectors = Point[]
    weights = []    
    for field in fields
        s, w = field.fieldFunction(p, field)
        push!(vectors, s)
        push!(weights, w)
    end
    weighted_sum_s = sum(weights[i] * vectors[i] for i in eachindex(vectors))
    total_w = sum(weights) + 0.01
    s_avg = weighted_sum_s / total_w

    s_avg /= (W/2)
    s_avg *= 20
    return s_avg
end

function attractorPoint(p::Point,F::attPointField)
    #rotation=0:no spiral, 1: ccw full rotation,-1 cw rotation
    #areaFactor: distance away at which weight is halfed
    p=p-F.c
    v=rotate(p,F.rotation*π/2)
    v*=.1
    weight=2^(-mag(v)^2/(F.areaFactor+.01))
    return v,weight
end



function attractorLine(p::Point, F::attLineField)
    AB = F.B - F.A
    Ap = p - F.A

    # Project p onto AB
    t = dotproduct(Ap, AB) / dotproduct(AB, AB)
    t = clamp(t, 0.0, 1.0)  # Clamp to segment
    closestPoint = F.A + t * AB

    d = distance(p, closestPoint)
    weight=2^(-d/(F.areaFactor+.01))

    gradient = (p - closestPoint)
    gradient=rotate(gradient,F.rotation*π/2)

    return gradient,weight
end

function lineSeeds(A::Point, B::Point, L, N)
    # Compute the direction of AB
    AB = B - A
    AB_norm = normalize(AB)

    n = Point(-AB_norm.y, AB_norm.x)

    # List of equally spaced points on the segment
    segment_points = [A + t * AB for t in range(0.0, stop=1.0, length=N÷2)]

    # Compute the offset points for both sides
    seeds = [point + L * n for point in segment_points]
    append!(seeds, [point - L * n for point in segment_points])
    # append!(seeds, [A - L * n, B + L * n])

    return seeds
end

function pointSeeds(c::Point, L, N)
    # List of equally spaced points in a polygon of radius L and N sides
    seeds = [Point(c.x+ L * cos(2 * π * i / N), c.y+ L * sin(2 * π * i / N)) for i in 1:N]

    return seeds
end


    
