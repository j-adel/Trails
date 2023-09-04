include("fns.jl")
include("Trail.jl")
include("fields.jl")
using Luxor

ntrails=100
W=H=600
trails=Vector{Trail}(undef,ntrails)



function main()
    trails = [Trail() for _=1:ntrails]  # Initialize an array of `ntrails` Trail objects with empty points arrays
    for trail in trails
        p=Point(rand(-W/2:W/2),rand(-H/2:H/2))
        followTrailBothWays!(trail,p)
        trail.origin=p
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