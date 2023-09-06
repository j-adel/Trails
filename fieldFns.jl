include("fns.jl")

function attractorPoint(p::Point,c::Point,rotation=0, areaFactor=Inf)
    #rotation=0:no spiral, 1: ccw full rotation,-1 cw rotation
    #areaFactor: distance away at which weight is halfed
    p=p-c
    v=rotate(p,rotation*π/2)
    v*=.1
    weight=2^(-mag(v)^2/(areaFactor+.01))
    return v,weight
end



function attractorLine(p::Point, A::Point, B::Point,
    rotation=0,areaFactor=Inf)
    AB = B - A
    Ap = p - A

    # Project p onto AB
    t = dotproduct(Ap, AB) / dotproduct(AB, AB)
    t = clamp(t, 0.0, 1.0)  # Clamp to segment
    closestPoint = A + t * AB

    d = distance(p, closestPoint)
    weight=2^(-d/(areaFactor+.01))

    gradient = (p - closestPoint)
    gradient=rotate(gradient,rotation*π/2)

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