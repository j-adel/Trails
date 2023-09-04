include("fns.jl")

function vectorfield(p)
    s0,w0=attractorLine(p,Point(0,-W/4),Point(0,W/4),.7,30)
    s1,w1=attractorLine(p,Point(-W/4,0),Point(W/4,0),.7,30)
    s=(w0*s0+w1+s1)/(w0+w1+.1)
    s/=(W/2)
    return s
end

function attractorPoint(p::Point,c::Point,rotate=0, areaFactor=Inf)
    #rotate=0:no spiral, 1: ccw full rotation,-1 cw rotation
    #areaFactor: distance away at which weight is halfed
    p=p-c
    v=rotate(p,rotate*π/2)
    v*=.1
    weight=2^(-mag(v)/(areaFactor+.01))
    return v,weight
end



function attractorLine(p::Point, A::Point, B::Point,
    rotationFactor=0,areaFactor=Inf)
    AB = B - A
    Ap = p - A

    # Project p onto AB
    t = dotproduct(Ap, AB) / dotproduct(AB, AB)
    t = clamp(t, 0.0, 1.0)  # Clamp to segment
    closestPoint = A + t * AB

    d = distance(p, closestPoint)
    weight=2^(-d/(areaFactor+.01))

    gradient = (p - closestPoint)
    gradient=rotate(gradient,rotationFactor*π/2)

    return gradient,weight
end




    
