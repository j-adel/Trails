using Luxor
include("fns.jl")

mutable struct Trail
    points::Vector{Point}
    origin::Point
    Trail() = new(Point[])
end

function followTrail!(trail,p,direction)
    maxNPts=2000
    i=0
    continueCond=true
    while p.x>-W/2 && p.x<W/2 && p.y>-H/2 && p.y<H/2 && continueCond
        s=vectorfield(p)
        p +=s*direction   # Update p1
        if direction==1
            push!(trail.points, p)  # Add the current p1 to t99rail.points
        elseif direction==-1
            pushfirst!(trail.points, p)  # Add the current p1 to start of trail.points
        else 
            error("direction must be 1 or -1")
        end
        if mag(s)<5/W && length(trail.points)>2
            if(direction==1)
                fp=findFixedPoint(trail.points[end-2],trail.points[end-1],trail.points[end])
                if !isnan(fp)
                    push!(trail.points,fp)
                    continueCond=false
                    println("found fixed point at ",fp," direction ",direction," after ",i," iterations")
                end
            else
                fp=findFixedPoint(trail.points[3],trail.points[2],trail.points[1])
                if !isnan(fp)
                    pushfirst!(trail.points,fp)
                    continueCond=false
                    println("found fixed point at ",fp," direction ",direction," after ",i," iterations")
                    println("trail.points=",trail.points[3],trail.points[2],trail.points[1])
                end
            end
        end
        i+=1
        continueCond=continueCond && i<maxNPts && mag(s)>1/W
    end
end

function findFixedPoint(A::Point,B::Point,C::Point)
    ratio=(C-B)/(B-A)
    # println("ratio=",ratio)
    # compute infinite geometric series
    return abs(ratio)<1 ? A + (B-A)/(1-ratio) : NaN
end



function followTrailBothWays!(trail,origin)
    p=origin #immutable
    push!(trail.points, p)
    followTrail!(trail,p,1)
    followTrail!(trail,p,-1)
end



function disp(trail)
    points=trail.points
    for i=1:length(points)-1
        line(points[i],points[i+1], :stroke)
    end
    gsave()
    sethue("red")
    circle(trail.origin,2,:fill)
    grestore()
end