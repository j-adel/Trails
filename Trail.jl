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
    maxNPts=1000
    i=0
    mergePulls=0
    contCond=true
    # p=Point(p.x,p.y+.5fieldFunction(p,fields).y)
    while contCond
        # RK2 integration
        s, continueCondSteps=integrateRK2(p,fieldFunction,direction,fields)
        contCond=contCond && continueCondSteps
        p += s
        #Verlet integration
        # p=Point(p.x+direction*fieldFunction(p,fields).x,p.y)
        # p=Point(p.x,p.y+direction*fieldFunction(p,fields).y)
        xIndex=limit(Int(round(p.x/cellSize)+ W÷(2*cellSize))+1,1,size(cells,1))
        yIndex=limit(Int(round(p.y/cellSize)+ H÷(2*cellSize))+1,1,size(cells,2))
        
        pNearest, dNearest=nearestPoint(p,trail,cells,xIndex,yIndex)
        # merge if close enough
        if dNearest< mergeDistance
            if dNearest<minStep
                p=pNearest
                contCond=false
                println("merged after ",i," iterations", " at ",p)
            else
                scaledAttractionCoeff=(.3+mergePulls/500)
                p=wavg(p,pNearest,max(0,scaledAttractionCoeff*parabola(dNearest,0,mergeDistance)))
                # println("moved ",p," closer to ",pNearest, "by", max(0,.3*parabola(dNearest,0,mergeDistance)))
                mergePulls+=1
            end
        else
            push!(cells[xIndex,yIndex],PointObject(p,trail))
        end

        # mergePulls>1 && println("Merge pulls: ",mergePulls)
        i+=1
        contCond=contCond && (length(points)<25|| mag(s)>minVectorLength || mag(s)>=mag(fieldFunction(points[end],fields)))
        if length(points)>25 && mag(p-points[end-1])/(mag(p-points[end])+mag(points[end]-points[end-1]))<.25
            # println("p ",p," points[end]",points[end]," points[end-1]",points[end-1])
            contCond=false    end
        if i>maxNPts-3 
            println("i: ", i, "p: ", p, " s: ",s) 
            push!(diagnosticPoints,p)
        end
        contCond=contCond && i<maxNPts && p.x>-.99W/2 && p.x<.99W/2 && p.y>-.99H/2 && p.y<.99H/2
        push!(points, p)  # Add the current p1 to trail.points
        # ~contCond && println("stopped after ",i," iterations", " at ",p, " with s ",s, " with continue conditions: ", 
        # " mag(s)>.05 ", mag(s)>.05, " mag(s)>=mag(fieldFunction(points[end],fields)) ", mag(s)>=mag(fieldFunction(points[end],fields)),
        # " mag ratio: ", length(points)>5 && mag(p-points[end])/(mag(p-points[end])+mag(points[end]-points[end-1]))<.25,
        # " i<maxNPts: ", i<maxNPts, " border collision: ", p.x>-.99W/2 && p.x<.99W/2 && p.y>-.99H/2 && p.y<.99H/2,"\n")

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

function integrateRK2(p::Point,fieldFunction::Function,direction,fields)
    s=Point(0,0)
    sMag=0
    steps=0
    maxSteps=minStep/minVectorLength*10
    while steps<maxSteps && mag(s)<minStep
        s+=fieldFunction(p+.5*direction*fieldFunction(p,fields),fields)*direction
        sMag+=mag(s)
        p += s
        steps+=1
    end
    steps>maxSteps-1 && println("steps: ",steps, " s: ",s)
    return s, mag(s)/sMag>.01
end

function nearestPoint(p,trail,cells,xIndex,yIndex)
    # find nearest neighbors
    dNearest=Inf
    pNearest=Point(Inf,Inf)
    for i=-1:1
        for j=-1:1
            xNbr=limit(xIndex+i,1,size(cells,1))
            yNbr=limit(yIndex+j,1,size(cells,2))
            for cellPointObject in cells[xNbr,yNbr]
                if cellPointObject.trail!=trail && cellPointObject.trail.sourceField!=trail.sourceField
                    dNearestTemp=distance(p,cellPointObject.point)
                    if dNearestTemp<dNearest
                        dNearest=dNearestTemp
                        pNearest=cellPointObject.point
                    end
                end
            end
        end
    end
    return pNearest, dNearest
end