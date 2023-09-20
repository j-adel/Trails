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
    areaFactor = 100 + 0abs(randn() * 50)
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