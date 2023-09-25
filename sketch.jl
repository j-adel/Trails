using Revise
using Plots
using Random
using Luxor
include("fns.jl")
include("Trail.jl")
include("Fields.jl")

W=H=600
cellSize=4 # integer
mergeDistance=cellSize
darkmode=false
displayColor="black"
trails=[]
distanceWPower=1.5; distancePower=.5
function main()
    darkmode ? displayColor="white" : displayColor="black"
    println(displayColor)
    SEED=rand(1:10000)
    SEED=4078
    Random.seed!(SEED)
    println("new run with seed: ",SEED)
    @time begin
        cells=Matrix{PointObject}(undef, W÷cellSize+2, H÷cellSize+2)
        fields=Field[]
        fieldsAmounts = Dict(attLineField => 2,streamLineField => 0, attPointField => 0, dipoleStreamField => 0)
        initializeFields!(fields,fieldsAmounts)
        println.(fields)
        _, _,trails=getSeeds(fields)
        for (i, trail) in enumerate(trails)
            # p = seeds[i]
            followTrailBothWays!(trail,computeVector,fields,cells)
            println("trail",i," length: ",length(trail.points))
        end
    end 
    

    @draw begin
        darkmode ? background("black") : background("white")
        for trail in trails
            disp(trail)
        end
        # for field in fields
        #     dispField(field)
        # end
        # #display seeds
        # sethue("red")
        # for trail in trails
        #     circle(trail.origin, 2, :fill)
        # end
        #display border
        sethue(displayColor)
        setline(2)
        rect(-W/2,-H/2, W, H, :stroke)
        Luxor.text(string(SEED), Point(-W/2-10,-H/2-10), halign=:left, valign=:bottom)
    end W+50 H+50
    # grid=range(-W/2,stop=W/2,length=100)
    # quiver_plot(computeVector,fields,20,20)
end

main()