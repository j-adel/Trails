using Luxor
include("fns.jl")
include("fieldFns.jl")

mutable struct Trail
    points::Vector{Point}
    origin::Point
    sourceField::Field
end

struct PointObject
    point::Point
    trail::Trail
end

function followTrail!(trail,p,direction,fieldFunction::Function,fields,cells)
    points=Point[]
    maxNPts=2000
    i=0
    contCond=true
    while p.x>-W/2 && p.x<W/2 && p.y>-H/2 && p.y<H/2 && contCond
        s=fieldFunction(p,fields)
        p +=s*direction   # Update p1
        xIndex=Int(round(p.x/cellSize)+ W÷(2*cellSize))+1; yIndex=Int(round(p.y/cellSize)+ H÷(2*cellSize))+1;
        # println(xIndex," ",yIndex)
        if isassigned(cells,xIndex,yIndex)
            if cells[xIndex,yIndex].trail!=trail && cells[xIndex,yIndex].trail.sourceField!=trail.sourceField
                d=distance(p,cells[xIndex,yIndex].point)
                if d<.5
                    p=cells[xIndex,yIndex].point
                    contCond=false
                else
                    p=wavg(p,cells[xIndex,yIndex].point,max(0,.5*parabola(d,0,cellSize)))
                end
            end
        else
            cells[xIndex,yIndex]=PointObject(p,trail)
        end
        push!(points, p)  # Add the current p1 to trail.points
        # if mag(s)<5/W && length(points)>5
        #     fp=findFixedPoint(points[end-2],points[end-1],points[end])
        #     if !isnan(fp)
        #         push!(points,fp)
        #         contCond=false
        #         # println("found fixed point at ",fp," direction ",direction," after ",i," iterations")
        #     end
        # end
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



function followTrailBothWays!(trail,origin,fieldFunction::Function,fields,cells)
    p=origin #immutable
    points1=followTrail!(trail,p,1,fieldFunction,fields,cells)
    points2=followTrail!(trail,p,-1,fieldFunction,fields,cells)
    append!(trail.points,reverse(points2))
    push!(trail.points, p)
    append!(trail.points,points1)
end



function disp(trail)
    if length(trail.points)<10
        return
    end
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