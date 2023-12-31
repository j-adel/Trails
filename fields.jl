
include("fns.jl")
include("fieldFns.jl")


function initializeFields!(fields::Vector{Field}, fieldsAmounts::Dict)
    for (fieldType, amount) in fieldsAmounts
        for _ = 1:amount
            push!(fields, initializeField(fieldType))
        end
    end
end



function getSeeds(fields,fieldFunction::Function)
    local seeds = Point[]
    seedTrails = Trail[]
    sourceFields = Field[]
    for field in fields
        if field.nSeeds > 0
            newSeedTrails = seedsFunction(field)(field,fieldFunction,fields)
            # append!(seeds, newSeeds)
            # append!(sourceFields, fill(field, length(newSeeds)))
            append!(seedTrails, newSeedTrails)
        end
    end
    return seedTrails
end

function getSeedsRnd(fields, nSeeds=50)
    seeds = Point[]
    sourceFields = Field[]
    for i = 1:nSeeds
        field = fields[rand(1:end)]
        newSeed = Point(randF(-W / 2, W / 2), randF(-H / 2, H / 2)) * 0.9
        push!(seeds, newSeed)
        append!(sourceFields, fill(field, length(newSeed)))
    end
    return seeds, sourceFields
end

function computeVector(p, fields)
    vectors = Point[]
    distances = []
    for field in fields
        s, w = fieldFunction(field)(field, p)
        push!(vectors, s)
        push!(distances, w)
    end
    # Inverse distance weighting
    sum_reciprocal = sum(1 / d^distanceWPower for d in distances)
    sum_vector_reciprocal = sum(vectors[i] / (distances[i]^distanceWPower) for i in eachindex(vectors))
    s_avg = sum_vector_reciprocal / sum_reciprocal

    s_avg *= 1
    if (abs(s_avg) > 30)
        println("s_avg: ", s_avg)
    end
    s_avg = containRect(p, s_avg, Point(-W / 2, -H / 2), Point(W / 2, H / 2))
    # mag(s_avg) > 10 && println("s_avg: ", s_avg)
    return s_avg
end


function quiver_plot(computeVector, fields, x_points::Int, y_points::Int)
    # Generate grid of points
    x_range, y_range = meshgrid(x_points, y_points, (-W / 2, W / 2), (-H / 2, H / 2))
    u = zeros(y_points, x_points)  # to store x-component of the vector
    v = zeros(y_points, x_points)  # to store y-component of the vector
    for i in 1:x_points
        for j in 1:y_points
            p = Point(x_range[i], y_range[j])
            vec = computeVector(p, fields)
            u[j, i] = vec[1]
            v[j, i] = vec[2]
        end
    end
    #print shape of u and v to make sure they are the same
    print(size(u))
    p = plot(; legend=false, showaxis=false)  # Creates an empty plot without legend and axis
    quiver!(p, x_range, y_range, quiver=(u, v))

end