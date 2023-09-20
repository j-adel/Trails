using Revise
using Plots
using Random
using Luxor
include("fns.jl")
include("Trail.jl")
include("Fields.jl")

W=H=600
cellSize=4 # integer
darkmode=false
trails=[]

function main()
    SEED=rand(1:10000)
    Random.seed!(SEED)
    @time begin
        cells=Matrix{PointObject}(undef, W÷cellSize+2, H÷cellSize+2)
        fields=Field[]
        fieldsAmounts = Dict(attLineField => 0,streamLineField => 0, attPointField => 0, dipoleStreamField => 1)
        initializeFields!(fields,fieldsAmounts)
        println.(fields)
        seeds, sourceFields=getSeedsRnd(fields)
        trails = [Trail(Point[],seeds[i],sourceFields[i]) for i in eachindex(seeds)]  # Initialize an  array of `ntrails` Trail objects with empty points arrays
        for (i, trail) in enumerate(trails)
            p = seeds[i]
            p=Point(rand(-W/2:W/2),rand(-H/2:H/2))
            followTrailBothWays!(trail, p,computeVector,fields,cells)
            trail.origin = p
        end
    end 
    

    @draw begin
        darkmode ? background("black") : background("white")
        for trail in trails
            disp(trail)
        end
        for field in fields
            dispField(field)
        end
        # draw a rectangle around the canvas
        sethue("black")
        setline(2)
        rect(-W/2,-H/2, W, H, :stroke)
        Luxor.text(string(SEED), Point(-W/2-10,-H/2-10), halign=:left, valign=:bottom)
    end W+50 H+50
    # # Display the SEED
    # pt1 = Point(-100, 0)
    # sethue("black")
    # quiver plot using Plots.jl of the vector field
    # grid=range(-W/2,stop=W/2,length=100)
    # quiver_plot(computeVector,fields,20,20)
end

main()