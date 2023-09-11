abstract type Field end
abstract type lineField <: Field end
abstract type pointField <: Field end

fieldFunction(f::Field) = error("fieldFunction not defined for type $f")
seedsFunction!(f::Field) = error("seedFunction not defined for type $f")

mutable struct attPointField <: pointField
    c::Point
    rotation::Float64
    areaFactor::Float64
    nSeeds::Int
    seedPoints::Vector{Point}
end
fieldFunction(::attPointField)=attPointFn
seedsFunction!(::attPointField)=pointSeedsFn!


function attPointFn(F::attPointField,p::Point)
    #rotation=0:no spiral, 1: ccw full rotation,-1 cw rotation
    #areaFactor: distance away at which weight is halfed
    v=p-F.c
    v=rotate(v,F.rotation*π/2)
    weight=2^(-mag(v)^2/(F.areaFactor+.01))
    v*=weight
    d=distance(p,F.c)/weight
    return v,d
end


function pointSeedsFn!(F::attPointField; N=10, L=0.1)
    seeds = Point[]
    for i in 1:N
        angle = rand() * 2π
        push!(seeds, F.c + Point(cos(angle), sin(angle)) * L)
    end
    return seeds
end



mutable struct attLineField <: lineField
    A::Point
    B::Point
    rotation::Float64
    areaFactor::Float64
    nSeeds::Int
    seedPoints::Vector{Point}
end
fieldFunction(::attLineField)=attLineFn
seedsFunction!(::attLineField)=lineSeedsFn!



function attLineFn(F::attLineField,p::Point)
    AB = F.B - F.A
    AB_norm = normalize(AB)
    Ap = p - F.A

    # Project p onto AB
    t = dotproduct(Ap, AB) / dotproduct(AB, AB)
    t = clamp(t, 0.0, 1.0)  # Clamp to segment
    closestPoint = F.A + t * AB

    gradient = (p - closestPoint)
    d = distance(p, closestPoint)
    weight=2^(-(d^2)/(F.areaFactor^2+.01))
    d /= weight
    gradient *= weight
    gradient=rotate(gradient,F.rotation*π/2)

    return gradient,d
end

mutable struct streamLineField <: lineField
    A::Point
    B::Point
    streamFactor::Float64
    areaFactor::Float64
    nSeeds::Int
    seedPoints::Vector{Point}
end
fieldFunction(::streamLineField)=streamLineFn
seedsFunction!(::streamLineField)=lineSeedsFn!

function streamLineFn(F::streamLineField,p::Point)
    AB_norm = normalize(F.B - F.A)

    # Distance to foci A and B
    dA = distance(p, F.A)
    dB = distance(p, F.B)

    # Calculate signed distance
    d = 0.5 * (dA + dB)

    # Normalize direction to foci A and B
    dirA = normalize(p - F.A)
    dirB = normalize(p - F.B)

    # The gradient is an average of the directions to the foci.
    # The 0.5 factor is because of the sum of distances.
    gradient = 0.5 * (dirA + dirB) * d
    weight=2^(-(d^2)/(F.areaFactor^2+.01))
    d /= weight
    gradient *= weight
    #stream factor weighting
    gradient=wavg(gradient,AB_norm*mag(gradient),map(weight,0,1,.8F.streamFactor,.9))

    return gradient,d
end

function lineSeedsFn!(F::lineField)
    # Compute the direction of AB
    L=0.1
    AB = F.B - F.A
    AB_norm = normalize(AB)

    n = Point(-AB_norm.y, AB_norm.x)
    segment_points = [F.A + t * AB for t in range(0.05, stop=.95, length=F.nSeeds÷2)]
    seeds = [point + L * n for point in segment_points]
    append!(seeds, [point - L * n for point in segment_points])
    # append!(seeds, [A - L * n, B + L * n])

    return seeds
end

function containRect(p::Point,v::Point,p1::Point,p2::Point;range=10)
    #shortest distance between p and the rectangle defined by p1 and p2
    
    d1=abs(p-p1)
    d2=abs(p-p2)
    d=min(d1,d2)
    if d>range
        return v
    end
    if d.x<d.y
        return Point(v.x*d.x/range,v.y*(d.x/range)^.5)
    else
        return Point(v.x*(d.y/range)^.5,v.y*d.y/range)
    end

end
