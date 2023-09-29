abstract type Field end

include("lineFieldFns.jl")
include("pointFieldFns.jl")
include("dipoleFieldFns.jl")

fieldFunction(f::Field) = error("fieldFunction not defined for type $f")
seedsFunction(f::Field) = error("seedFunction not defined for type $f")

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






