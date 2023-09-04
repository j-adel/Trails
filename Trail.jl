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
        i+=1
        continueCond=i<maxNPts && mag(s)>1/W
    end
end

function followTrailBothWays!(trail,origin=Point(rand(-w/2:w/2),rand(-h/2:h/2)))
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
    # #draw point at trail.origin in red
    # gsave()
    # sethue("red")
    # circle(trail.origin,2,:fill)
    # grestore()
end