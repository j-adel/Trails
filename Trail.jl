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


function followTrail(trail,p,direction,fieldFunction::Function,fields,cells)
    points=Point[]
    maxNPts=5000
    i=0
    contCond=true
    # p=Point(p.x,p.y+.5fieldFunction(p,fields).y)
    while contCond
        s=fieldFunction(p,fields)*direction
        # RK2 integration
        p +=fieldFunction(p+.5s,fields)*direction   # Update p1
        #Verlet integration
        # p=Point(p.x+direction*fieldFunction(p,fields).x,p.y)
        # p=Point(p.x,p.y+direction*fieldFunction(p,fields).y)
        xIndex=limit(Int(round(p.x/cellSize)+ W÷(2*cellSize))+1,1,size(cells,1))
        yIndex=limit(Int(round(p.y/cellSize)+ H÷(2*cellSize))+1,1,size(cells,2))
        
        # find nearest neighbors
        dNearest=Inf
        pNearest=Point(Inf,Inf)
        for i=-1:1
            for j=-1:1
                xNbr=limit(xIndex+i,1,size(cells,1))
                yNbr=limit(yIndex+j,1,size(cells,2))
                if isassigned(cells,xNbr,yNbr)
                    if cells[xNbr,yNbr].trail!=trail #&& cells[xNbr,yNbr].trail.sourceField!=trail.sourceField
                        dNearestTemp=distance(p,cells[xNbr,yNbr].point)
                        if dNearestTemp<dNearest
                            dNearest=dNearestTemp
                            pNearest=cells[xNbr,yNbr].point
                        end
                    end
                end
            end
        end
        # merge if close enough
        if dNearest< mergeDistance
            if dNearest<.5
                p=pNearest
                contCond=false
                # println("merged after ",i," iterations", " at ",p)
            else
                p=wavg(p,pNearest,max(0,.3*parabola(dNearest,0,mergeDistance)))
            end
        elseif ~isassigned(cells,xIndex,yIndex)
            cells[xIndex,yIndex]=PointObject(p,trail)
        end
        push!(points, p)  # Add the current p1 to trail.points
        i+=1
        contCond=contCond && (mag(s)>.1/W || mag(s)>mag(fieldFunction(p-s,fields)) || length(points)<5)
        contCond=contCond && i<maxNPts && p.x>-W/2 && p.x<W/2 && p.y>-H/2 && p.y<H/2
        # ~contCond && println("stopped after ",i," iterations", " at ",p, " with s ",s)
    end
    return points
end



function followTrailBothWays!(trail,fieldFunction::Function,fields,cells)
    #find the direction of the Trail
    direction1=sign(dot(fieldFunction(trail.points[2],fields),trail.points[2]-trail.origin))
    direction2=sign(dot(fieldFunction(trail.points[1],fields),trail.points[1]-trail.origin))
    points1=followTrail(trail,trail.points[end],direction1,fieldFunction,fields,cells)
    points2=followTrail(trail,trail.points[1],direction2,fieldFunction,fields,cells)
    trail.points=append!(reverse(points2),trail.points)
    append!(trail.points,points1)
end


function findFixedPoint(A::Point,B::Point,C::Point)
    ratio=(C-B)/(B-A)
    # compute infinite geometric series
    return abs(ratio)<.999 ? A + (B-A)/(1-ratio) : NaN
end


function disp(trail)
    if length(trail.points)<10
        return
    end
    points=trail.points
    move(points[1])
    setcolor(displayColor)
    setline(1.5)
    for i=1+1:length(points)
        line(points[i])
        # gsave()
        # sethue("red")
        # circle(points[i],2,:fill)
        # grestore()
    end
    strokepath()
    
    # gsave()
    # sethue("red")
    # circle(trail.origin,2,:fill)
    # grestore()
end

# function integrate(p::Point,fieldFunction::Function,fields)
#     p=Point(fieldFunction(p,fields).x,p.y)
#     p=Point(p.x,fieldFunction(p,fields).y)
#     return p
# end