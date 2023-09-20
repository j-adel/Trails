################ dipole Line Field ################
abstract type dipoleField <: Field end
mutable struct dipoleStreamField <: dipoleField
    P::Point
    direction::Float64
    dipoleFactor::Float64
    areaFactor::Float64
    nSeeds::Int
    seedPoints::Vector{Point}
end
fieldFunction(::dipoleStreamField)=dipoleStreamFn
seedsFunction!(::dipoleStreamField)=error("dipoleStreamField does not have seeds")

function initializeField(::Type{T}) where T <: dipoleStreamField
    P=Point(randF(-W/2, W/2), randF(-H/2, H/2))*.8
    direction = rand()*2π
    dipoleFactor = rand()
    areaFactor = 1
    nSeeds = 0
    seedPoints = Point[]
    
    return T(P, direction, dipoleFactor, areaFactor, nSeeds, seedPoints)
end

function dipoleStreamFn(F::dipoleStreamField,p::Point)
    d=distance(p,F.P)
    directionVec=Point(cos(F.direction),sin(F.direction))
    # convexity=dot(p-F.P,directionVec)/(1W) #radius at d=1
    # # println("convexity: ",convexity)
    # c=F.P+directionVec/sP(convexity,3)
    # gradient=perpendicular(norm(p-c))/(d+1)*W/2
    # gradient=rotate((p-F.P),-smoothStep(heading(p-F.P)/π)*π)
    angle=modulo(heading(p-F.P)+F.direction,2π)-π
    if angle>0
        gradient=rotate((p-F.P),-scaledStep(angle,-d/W*5,0,π,0,π))
    else
        gradient=rotate((p-F.P),scaledStep(-angle,-d/W*5,0,π,0,π))
    end

    return gradient,d
end

mutable struct dipoleRingField <: dipoleField
    A::Point
    B::Point
    dipoleFactor::Float64
    areaFactor::Float64
    nSeeds::Int
    seedPoints::Vector{Point}
end