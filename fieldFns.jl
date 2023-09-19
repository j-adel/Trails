abstract type Field end

fieldFunction(f::Field) = error("fieldFunction not defined for type $f")
seedsFunction!(f::Field) = error("seedFunction not defined for type $f")

################ Attractive Point Field ################
abstract type pointField <: Field end
mutable struct attPointField <: pointField
    c::Point
    rotation::Float64
    areaFactor::Float64
    nSeeds::Int
    seedPoints::Vector{Point}
end
fieldFunction(::attPointField)=attPointFn
seedsFunction!(::attPointField)=pointSeedsFn!

function initializeField(::Type{T}) where T <: attPointField
    center = Point(randF(-W/2, W/2), randF(-H/2, H/2))
    rotation = randF(-1, 1)*.9 +rand([0,2])
    areaFactor = abs(randn() * 50)+50
    nSeeds = randI(5,12)
    seedPoints = Point[]
    
    return T(center, rotation, areaFactor, nSeeds, seedPoints)
end

function attPointFn(F::attPointField,p::Point)
    #rotation=0:no spiral, 1: ccw full rotation,-1 cw rotation
    #areaFactor: distance away at which weight is halfed
    v=p-F.c
    v=rotate(v,F.rotation*π/2)
    weight=2^(-mag(v)^2/(F.areaFactor^2+.01))
    v*=weight
    d=distance(p,F.c)/weight
    return v,d
end

function pointSeedsFn!(F::attPointField; )
    seeds = Point[]
    for i in 1:F.nSeeds
        angle = map(i, 0, F.nSeeds, 0, 2π)
        push!(seeds, F.c + Point(cos(angle), sin(angle)) * 1)
    end
    return seeds
end

################ Attractive Line Field ################

abstract type lineField <: Field end
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

function initializeField(::Type{T}) where T <: attLineField
    center = Point(randF(-W/3, W/3), randF(-H/3, H/3))
    angle = rand() * 2π
    length = randF(0.1, 0.5) * W
    
    A = limitMag(center + Point(cos(angle), sin(angle)) * length, Point(W/2, H/2)*.9)
    B = limitMag(center + Point(cos(angle + π), sin(angle + π)) * length, Point(W/2, H/2)*.9)
    rotation = randF(-1, 1)*.9 +rand([-2,0,2])
    areaFactor = abs(randn() * 50)+distance(A,B)/4
    nSeeds = 14
    seedPoints = Point[]
    
    return T(A, B, rotation, areaFactor, nSeeds, seedPoints)
end

function attLineFn(F::attLineField,p::Point)
    AB = F.B - F.A
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

################ Stream Line Field ################

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

function initializeField(::Type{T}) where T <: streamLineField
    center = Point(randF(-W/3, W/3), randF(-H/3, H/3))
    angle = rand() * 2π
    length = randF(0.2, 0.5) * W
    
    A = limitMag(center + Point(cos(angle), sin(angle)) * length, Point(W/2, H/2)*.9)
    B = limitMag(center + Point(cos(angle + π), sin(angle + π)) * length, Point(-W/2, H/2)*.9)
    streamFactor = rand()
    areaFactor = abs(randn() * 50)+distance(A,B)/3
    nSeeds = 0
    seedPoints = Point[]
    
    return T(A, B, streamFactor, areaFactor, nSeeds, seedPoints)
end

function streamLineFn(F::streamLineField,p::Point)
    AB_norm = normalize(F.B - F.A)

    gradient,d = ellipseGrad(p, F.A, F.B)
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
    segment_points = [F.A + t * AB for t in range(0.0, stop=1, length=F.nSeeds÷2)]
    seeds = [point + L * n for point in segment_points]
    append!(seeds, [point - L * n for point in segment_points])
    # append!(seeds, [A - L * n, B + L * n])

    return seeds
end

function containRect(p::Point,v::Point,p1::Point,p2::Point;range=20)
    #shortest distance between p and the rectangle defined by p1 and p2
    
    d1=abs(p-p1)
    d2=abs(p-p2)
    d=min(d1,d2)
    if d>range
        return v
    end
    if d.x<d.y
        return Point(v.x*d.x/range,v.y*(d.x/range)^.25)
    else
        return Point(v.x*(d.y/range)^.25,v.y*d.y/range)
    end

end

function ellipseGrad(p::Point, A::Point, B::Point)
    # Distance to foci A and B
    dA = distance(p, A)
    dB = distance(p, B)

    # Calculate signed distance
    d = 0.5 * (dA + dB)

    # Normalize direction to foci A and B
    dirA = normalize(p - A)
    dirB = normalize(p - B)

    # The gradient is an average of the directions to the foci.
    # The 0.5 factor is because of the sum of distances.
    gradient = 0.5 * (dirA + dirB) * d
    return gradient,d
end



