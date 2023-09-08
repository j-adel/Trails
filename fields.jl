include("fns.jl")
abstract type Field end

fieldFunction(f::Field) = error("fieldFunction not defined for type $f")
seedsFunction!(f::Field) = error("seedFunction not defined for type $f")

mutable struct attPointField <: Field
    c::Point
    rotation::Float64
    areaFactor::Float64
    nSeeds::Int
    seedPoints::Vector{Point}
    attPointField() = new(undef, undef, undef, undef, undef)
end
fieldFunction(::attPointField)=attPointFn
seedsFunction!(::attPointField)=pointSeedsFn!

mutable struct attLineField <: Field
    A::Point
    B::Point
    rotation::Float64
    areaFactor::Float64
    nSeeds::Int
    seedPoints::Vector{Point}
end
fieldFunction(::attLineField)=attLineFn
seedsFunction!(::attLineField)=lineSeedsFn!

function initializeFields!(fields::Vector{Field}, nfields=2)
    for i=1:nfields
        field=initializeField(attLineField)
        push!(fields, field)
    end

    # push!(fields,attLineField(Point(-W/2,-H/2),Point( W/2,-H/2),0.95,0,0,Point[]))
    # push!(fields,attLineField(Point( W/2,-H/2),Point( W/2, H/2),0.95,0,0,Point[]))
    # push!(fields,attLineField(Point( W/2, H/2),Point(-W/2, H/2),0.95,0,0,Point[]))
    # push!(fields,attLineField(Point(-W/2, H/2),Point(-W/2,-H/2),0.95,0,0,Point[]))

    # field1=attLineField(Point(0, -W/4), Point(0, W/4), 1.5, 50, 10, Point[])
    # push!(fields, field1)
    # field2=attLineField(Point(-W/4, 0), Point(W/4, 0), 0.3, 50, 10,Point[])
    # push!(fields, field2)
    # field3=attPointField(Point(0, 0), 2, 150, Point[])
    # push!(fields, field3)
    # push!(fields, attPointField((p)->vectorfield(p), Point(0, 0), 0.5, 150, Point[]))
end

function initializeField(::Type{T}) where T <: attLineField
    center = Point(randF(-W/3, W/3), randF(-H/3, H/3))
    angle = rand() * 2π
    length = randF(0.1, 0.5) * W
    
    A = max(min(center + Point(cos(angle), sin(angle)) * length, Point(W/2, H/2)*.9), Point(-W/2, -H/2)*.9)
    B = max(min(center + Point(cos(angle + π), sin(angle + π)) * length, Point(W/2, H/2)*.9), Point(-W/2, -H/2)*.9)
    rotation = randF(-1, 1)*.95 +rand([-1,0,1])
    areaFactor = abs(randn() * 50)+distance(A,B)/4
    nSeeds = 20
    seedPoints = Point[]
    
    return T(A, B, rotation, areaFactor, nSeeds, seedPoints)
end



function getSeeds(fields)
    seeds=Point[]
    for field in fields
        append!(seeds,seedsFunction!(field)(field))
    end
    return seeds
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
    if ~(abs(s_avg)<30)
        println("s_avg: ",s_avg)
    end
    return s_avg
end

function attPointFn(F::attPointField,p::Point)
    #rotation=0:no spiral, 1: ccw full rotation,-1 cw rotation
    #areaFactor: distance away at which weight is halfed
    v=p-F.c
    v=rotate(v,F.rotation*π/2)
    weight=2^(-mag(v)^2/(F.areaFactor+.01))
    v*=weight
    d=distance(p,F.c)/weight
    return v,d
end



function attLineFn(F::attLineField,p::Point)
    AB = F.B - F.A
    Ap = p - F.A

    # Project p onto AB
    t = dotproduct(Ap, AB) / dotproduct(AB, AB)
    t = clamp(t, 0.0, 1.0)  # Clamp to segment
    closestPoint = F.A + t * AB

    gradient = (p - closestPoint)
    d = distance(p, closestPoint)
    weight=2^(-d/(F.areaFactor+.01))
    d /= weight
    gradient *= weight
    gradient=rotate(gradient,F.rotation*π/2)

    return gradient,d
end

function lineSeedsFn!(F::attLineField; N=10, L=0.1)
    # Compute the direction of AB
    AB = F.B - F.A
    AB_norm = normalize(AB)

    n = Point(-AB_norm.y, AB_norm.x)
    segment_points = [F.A + t * AB for t in range(0.0, stop=1.0, length=N÷2+1)]
    seeds = [point + L * n for point in segment_points[1:end-1]]
    append!(seeds, [point - L * n for point in segment_points[2:end]])
    # append!(seeds, [A - L * n, B + L * n])

    return seeds
end

function pointSeedsFn!(F::attPointField; N=10, L=0.1)
    seeds = Point[]
    for i in 1:N
        angle = rand() * 2π
        push!(seeds, F.c + Point(cos(angle), sin(angle)) * L)
    end
    return seeds
end

function dispField(F::attPointField)
    # display the element
    gsave()
    darkmode ? sethue("white") : sethue("black")
    circle(F.c,2,:fill)
    grestore()
end

function dispField(F::attLineField)
    # display the element
    gsave()
    setline(1.5)
    darkmode ? sethue("white") : sethue("black")
    # sethue(0,0,.5)
    line(F.A,F.B,:stroke)
    grestore()
end
    
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