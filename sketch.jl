using Revise
using Plots
using Random
using Luxor
include("fns.jl")
include("Trail.jl")
include("Fields.jl")

W=H=600
cellSize=4 # integer
minVectorLength=.01; minStep=.5
mergeDistance=cellSize
darkmode=false
displayColor="black"
diagnosticMode=false
diagnosticPoints=Point[]
trails=[]
distanceWPower=1.5; distancePower=.5
function main()
    darkmode ? displayColor="white" : displayColor="black"
    SEED=rand(1:10000)
    # SEED=2650
    Random.seed!(SEED)
    println("new run with seed: ",SEED)
    @time begin
        cells = [PointObject[] for _ in 1:(W ÷ cellSize + 2), _ in 1:(H ÷ cellSize + 2)]
        fields=Field[]
        fieldsAmounts = Dict(attLineField =>3,streamLineField => 0, attPointField => 0, dipoleStreamField => 0)
        initializeFields!(fields,fieldsAmounts)
        println.(fields)
        trails=getSeeds(fields,computeVector)
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
        

        if diagnosticMode #display seeds
            sethue("red")
            # for trail in trails
            #     circle(trail.origin, 2, :fill)
            # end
            for field in fields
                dispField(field)
            end
            for point in diagnosticPoints
                circle(point, 2, :stroke)
            end
        end
        #display border
        sethue(displayColor)
        setline(3)
        rect(-.99W/2,-.99H/2, .99W, .99H, :stroke)
        Luxor.text(string(SEED), Point(-W/2-10,-H/2-10), halign=:left, valign=:bottom)
    end W+50 H+50
    # grid=range(-W/2,stop=W/2,length=100)
    # quiver_plot(computeVector,fields,20,20)
end

main()