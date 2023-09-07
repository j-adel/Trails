using Luxor
include("fns.jl")
include("Trail.jl")
include("fields.jl")

ntrails=300
nfields=5
W=H=600
trails=Vector{Trail}(undef,ntrails)



function main()
    @time begin
    trails = [Trail() for _=1:ntrails]  # Initialize an array of `ntrails` Trail objects with empty points arrays
    seeds=lineSeeds(Point(0, -W/4), Point(0, W/4),15,ntrails)
    fields=[]
    initializeFields!(fields)
    # append!(seeds,lineSeeds(Point(-W/4, 0), Point(W/4, 0),15,ntrails÷2))
    for (i, trail) in enumerate(trails)
        # p = seeds[i]
        p=Point(rand(-W/2:W/2),rand(-H/2:H/2))
        followTrailBothWays!(trail, p,computeVector,fields)
        trail.origin = p
    end
    

    @draw begin
        background("black")
        for trail in trails
            disp(trail)
        end
    end W H
    end 
end

main()