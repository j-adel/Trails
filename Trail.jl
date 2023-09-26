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
    mergePulls=0
    contCond=true
    # p=Point(p.x,p.y+.5fieldFunction(p,fields).y)
    while contCond
        # RK2 integration
        s=fieldFunction(p+.5*direction*fieldFunction(p,fields),fields)*direction
        p += s  # Update p1
        #Verlet integration
        # p=Point(p.x+direction*fieldFunction(p,fields).x,p.y)
        # p=Point(p.x,p.y+direction*fieldFunction(p,fields).y)
        xIndex=limit(Int(round(p.x/cellSize)+ W÷(2*cellSize))+1,1,size(cells,1))
        yIndex=limit(Int(round(p.y/cellSize)+ H÷(2*cellSize))+1,1,size(cells,2))
        
        # find nearest neighbors
        dNearest=Inf
        pNearest=Point(Inf,-Inf)
        for i=-1:1
            for j=-1:1
                xNbr=limit(xIndex+i,1,size(cells,1))
                yNbr=limit(yIndex+j,1,size(cells,2))
                if isassigned(cells,xNbr,yNbr)
                    if cells[xNbr,yNbr].trail!=trail && cells[xNbr,yNbr].trail.sourceField!=trail.sourceField
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
                println("merged ", p, " with ",pNearest)
                p=pNearest
                contCond=false
                # println("merged after ",i," iterations", " at ",p)
            else
                scaledAttractionCoeff=(.3+mergePulls/500)
                p=wavg(p,pNearest,max(0,scaledAttractionCoeff*parabola(dNearest,0,mergeDistance)))
                # println("moved ",p," closer to ",pNearest, "by", max(0,.3*parabola(dNearest,0,mergeDistance)))
                mergePulls+=1
            end
        elseif ~isassigned(cells,xIndex,yIndex)
            cells[xIndex,yIndex]=PointObject(p,trail)
        end

        # mergePulls>1 && println("Merge pulls: ",mergePulls)
        push!(points, p)  # Add the current p1 to trail.points
        i+=1
        contCond=contCond && (mag(s)>.01 || mag(s)>=mag(fieldFunction(p-s,fields)) || length(points)<25)
        contCond=contCond && i<maxNPts && p.x>-W/2 && p.x<W/2 && p.y>-H/2 && p.y<H/2
        # ~contCond && println("stopped after ",i," iterations", " at ",p, " with s ",s)
    end
    return points
end



function followTrailBothWays!(trail,fieldFunction::Function,fields,cells)
    #find the direction of the Trail
    lookout=0
    # println("points",trail.points[1],trail.points[2],wavg(trail.points[1],trail.points[2],1+lookout))
    direction1=sign(dot(fieldFunction(wavg(trail.points[1],trail.points[2],0-lookout),fields),trail.points[1]-trail.origin))
    direction2=sign(dot(fieldFunction(wavg(trail.points[1],trail.points[2],1+lookout),fields),trail.points[2]-trail.origin))
    # println("directions ",direction1,direction2)
    # direction1=1;direction2=1;
    points1=followTrail(trail,trail.points[1],direction1,fieldFunction,fields,cells)
    points2=followTrail(trail,trail.points[end],direction2,fieldFunction,fields,cells)
    trail.points=append!(reverse(points1),[trail.points[1],trail.points[2]])
    append!(trail.points,points2)
end


function findFixedPoint(A::Point,B::Point,C::Point)
    ratio=(C-B)/(B-A)
    # compute infinite geometric series
    return abs(ratio)<.999 ? A + (B-A)/(1-ratio) : NaN
end


function disp(trail)
    if length(trail.points)<10
        return
    end ######### CHANGE THIS TO DISTANCE
    points=trail.points
    setcolor(displayColor)
    setline(1.5)
    bezpath=makebezierpath(points)
    drawbezierpath(bezpath,:stroke,close=false)
    # move(points[1])
    # for i=1+1:length(points)
    #     line(points[i])
    #     # gsave()
    #     # sethue("red")
    #     # circle(points[i],2,:fill)
    #     # grestore()
    # end
    # strokepath()
    
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