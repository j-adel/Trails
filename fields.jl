
include("fns.jl")
include("fieldFns.jl")


function initializeFields!(fields::Vector{Field}, fieldsAmounts::Dict)
    for (fieldType, amount) in fieldsAmounts
        for _=1:amount
            push!(fields, initializeField(fieldType))
        end
    end
end



function getSeeds(fields)
    local seeds=Point[]
    sourceFields=Field[]
    for field in fields
        if field.nSeeds>0
            newSeeds=seedsFunction!(field)(field)
            append!(seeds,newSeeds)
            append!(sourceFields,fill(field,length(newSeeds)))
        end
    end
    return seeds, sourceFields
end

function getSeedsRnd(fields, nSeeds=50)
    seeds=Point[]
    sourceFields=Field[]
    for i=1:nSeeds
        field=fields[rand(1:end)]
        newSeed=Point(randF(-W/2,W/2),randF(-H/2,H/2))*.9
        push!(seeds,newSeed)
        append!(sourceFields,fill(field,length(newSeed)))
    end
    return seeds, sourceFields
end

function computeVector(p,fields)
    vectors = Point[]
    distances = []    
    for field in fields
        s, w = fieldFunction(field)(field, p)
        push!(vectors, s)
        push!(distances, w)
    end
    # Inverse distance weighting
    sum_reciprocal = sum(1 / d for d in distances)
    sum_vector_reciprocal = sum(vectors[i] / distances[i] for i in eachindex(vectors))
    s_avg = sum_vector_reciprocal / sum_reciprocal
    
    s_avg /= (W/2)
    s_avg *= 20
    if (abs(s_avg)>30)
        println("s_avg: ",s_avg)
    end
    s_avg=containRect(p,s_avg,Point(-W/2,-H/2),Point(W/2,H/2))
    return s_avg
end

function dispField(F::attPointField)
    if F.nSeeds==0
        return
    end
    gsave()
    darkmode ? sethue("white") : sethue("black")
    circle(F.c,2,:fill)
    grestore()
end

function dispField(F::lineField)
    # display the element
    gsave()
    setline(1.5)
    darkmode ? sethue("white") : sethue("black")
    # sethue(0,0,.5)
    line(F.A,F.B,:stroke)
    grestore()
end
#dispField of dipoleField is empty
dispField(F::dipoleField)=nothing
    
function quiver_plot(computeVector,fields, x_points::Int, y_points::Int)
    # Generate grid of points
    x_range, y_range = meshgrid(x_points, y_points,(-W/2,W/2), (-H/2,H/2))
    u = zeros(y_points, x_points)  # to store x-component of the vector
    v = zeros(y_points, x_points)  # to store y-component of the vector
    for i in 1:x_points
        for j in 1:y_points
            p = Point(x_range[i], y_range[j])
            vec = computeVector(p,fields)
            u[j, i] = vec[1]
            v[j, i] = vec[2]
        end
    end
    #print shape of u and v to make sure they are the same
    print(size(u))
    p = plot(; legend=false, showaxis=false)  # Creates an empty plot without legend and axis
    quiver!(p, x_range, y_range, quiver=(u, v))
    
end