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
    rotation = randF(-1, 1)*.8 +rand([-2,0,2])
    areaFactor = distance(A,B)/2 + 0abs(randn() * 50)
    nSeeds = 14
    seedPoints = Point[]
    
    return T(A, B, rotation, areaFactor, nSeeds, seedPoints)
end

function attLineFn(F::attLineField,p::Point)
    # AB = F.B - F.A
    # Ap = p - F.A

    # # Project p onto AB
    # t = dotproduct(Ap, AB) / dotproduct(AB, AB)
    # t = clamp(t, 0.0, 1.0)  # Clamp to segment
    # closestPoint = F.A + t * AB

    # gradient = (p - closestPoint)
    # d = distance(p, closestPoint)
    gradient,d = ellipseGrad(p, F.A, F.B)
    weight=2^(-(d^2)/(F.areaFactor^2+.01))
    # d /= weight
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
    areaFactor = distance(A,B)/2+ 0abs(randn() * 50)
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
    # gradient=wavg(gradient,AB_norm*mag(gradient),map(weight,0,1,.8F.streamFactor,.9))
    gradient=AB_norm*mag(gradient)
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