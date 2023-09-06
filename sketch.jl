using Luxor
include("fns.jl")
include("Trail.jl")
include("fields.jl")

ntrails=12
W=H=600
trails=Vector{Trail}(undef,ntrails)



function main()
    trails = [Trail() for _=1:ntrails]  # Initialize an array of `ntrails` Trail objects with empty points arrays
    seeds=lineSeeds(Point(0, -W/4), Point(0, W/4),15,ntrails)
    # append!(seeds,lineSeeds(Point(-W/4, 0), Point(W/4, 0),15,ntrails÷2))
    for (i, trail) in enumerate(trails)
        # p = seeds[i]
        p=Point(rand(-W/2:W/2),rand(-H/2:H/2))
        followTrailBothWays!(trail, p)
        trail.origin = p
    end
    

    @draw begin
        background("black")
        setcolor(1,1,1)
        for trail in trails
            disp(trail)
        end
    
    end W H
end

main()