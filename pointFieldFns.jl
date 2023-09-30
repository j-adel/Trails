################ Attractive Point Field ################

abstract type pointField <: Field end
mutable struct attPointField <: pointField
    c::Point
    rotation::Float64
    direction::Int8
    areaFactor::Float64
    nSeeds::Int
    seedPoints::Vector{Point}
end
fieldFunction(::attPointField)=attPointFn
seedsFunction(::attPointField)=pointSeedsFn

function initializeField(::Type{T}) where T <: attPointField
    center = Point(randF(-W/2, W/2), randF(-H/2, H/2))
    rotation = randF(-1, 1)*.8
    direction=rand([-1,1])
    areaFactor = 100 + 0abs(randn() * 50)
    nSeeds = randI(5,12)
    seedPoints = Point[]
    
    return T(center, rotation, direction, areaFactor, nSeeds, seedPoints)
end

function attPointFn(F::attPointField,p::Point)
    #rotation=0:no spiral, 1: ccw full rotation,-1 cw rotation
    #areaFactor: distance away at which weight is halfed
    v=p-F.c
    d=distance(p,F.c)
    v=rotate(v,F.rotation*π/2)*((d^distancePower)/d)*F.direction
    weight=2^(-mag(v)^2/(F.areaFactor^2+.01))
    # v*=weight
    d=distance(p,F.c)/weight
    return v,d
end

function pointSeedsFn(F::attPointField,fieldFunction::Function,fields::Vector{Field})
    L = 1 # You can adjust the length L to your needs
    seedTrails = Trail[]
    for i in 1:F.nSeeds
        angle = map(i, 1, F.nSeeds, 0, 2π)
        point1 = F.c + Point(cos(angle), sin(angle)) * L
        point2 = F.c + Point(cos(angle + π), sin(angle + π)) * L # point at the opposite end of the circle
        push!(seedTrails, Trail([point1, point2], F.c, F))
    end
    return seedTrails
end

function dispField(F::attPointField)
    if F.nSeeds == 0
        return
    end
    gsave()
    darkmode ? sethue("white") : sethue("black")
    circle(F.c, 2, :fill)
    grestore()
end

