using Luxor
using Plots
include("fns.jl")
include("Trail.jl")
include("Fields.jl")

W=H=600
darkmode=false
trails=[]


function main()
    @time begin
        # append!(seeds,lineSeedsFn(Point(-W/4, 0), Point(W/4, 0),1,ntrails÷2))
        fields=Field[]
        fieldsAmounts = Dict(attPointField => 0, attLineField => 0, streamLineField =>3)
        initializeFields!(fields,fieldsAmounts)
        println.(fields)
        seeds=getSeeds(fields)
        trails = [Trail() for _=1:length(seeds)]  # Initialize an  array of `ntrails` Trail objects with empty points arrays
        for (i, trail) in enumerate(trails)
            p = seeds[i]
            # p=Point(rand(-W/2:W/2),rand(-H/2:H/2))
            followTrailBothWays!(trail, p,computeVector,fields)
            trail.origin = p
        end
    

    @draw begin
        darkmode ? background("black") : background("white")
        for trail in trails
            disp(trail)
        end
        for field in fields
            dispField(field)
        end
    end W+50 H+50
    # quiver plot using Plots.jl of the vector field
    # grid=range(-W/2,stop=W/2,length=100)
    # quiver_plot(computeVector,fields,20,20)
    end 
end

main()