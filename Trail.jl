using Luxor
include("fns.jl")

mutable struct Trail
    points::Vector{Point}
    origin::Point
    Trail() = new(Point[])
end

function followTrail!(p,direction,fieldFunction::Function,fields)
    points=Point[]
    maxNPts=2000
    i=0
    contCond=true
    while p.x>-W/2 && p.x<W/2 && p.y>-H/2 && p.y<H/2 && contCond
        s=fieldFunction(p,fields)
        p +=s*direction   # Update p1
        push!(points, p)  # Add the current p1 to trail.points
        if mag(s)<5/W && length(points)>5
            fp=findFixedPoint(points[end-2],points[end-1],points[end])
            if !isnan(fp)
                push!(points,fp)
                contCond=false
                # println("found fixed point at ",fp," direction ",direction," after ",i," iterations")
            end
        end
        i+=1
        contCond=contCond && i<maxNPts && (mag(s)>.1/W || length(points)<5)
    end
    return points
end

function findFixedPoint(A::Point,B::Point,C::Point)
    ratio=(C-B)/(B-A)
    # compute infinite geometric series
    return abs(ratio)<.999 ? A + (B-A)/(1-ratio) : NaN
end



function followTrailBothWays!(trail,origin,fieldFunction::Function,fields)
    p=origin #immutable
    points1=followTrail!(p,1,fieldFunction,fields)
    points2=followTrail!(p,-1,fieldFunction,fields)
    append!(trail.points,reverse(points2))
    push!(trail.points, p)
    append!(trail.points,points1)
end



function disp(trail)
    points=trail.points
    move(points[1])
    darkmode ? setcolor("white") : setcolor("black")
    setline(1.5)
    for i=1+1:length(points)
        line(points[i])
    end
    strokepath()
    
    # gsave()
    # sethue("red")
    # circle(trail.origin,2,:fill)
    # grestore()
end