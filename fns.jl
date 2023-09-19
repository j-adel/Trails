using Luxor

function mag(p::Point)
    return sqrt(p.x^2+p.y^2)
end

function perpendicular(p::Point)
    return Point(-p.y,p.x)
end

function map(x,x0,x1,y0,y1)
    return (x-x0)/(x1-x0)*(y1-y0)+y0
end

function mapto01(x,x0,x1)
    return (x-x0)/(x1-x0)
end

function mapfrom01(x,y0,y1)
    return x*(y1-y0)+y0
end

function wavg(x0,x1,w=.5)
    return (1-w)*x0+w*x1
end


function step(x,x0=0)
    return x>x0 ? 1 : 0
end

function normalize(p)
    len = sqrt(p.x^2 + p.y^2)
    if len != 0
        p = Point(p.x / len, p.y / len)
    end
    return p
end

function rotate(p,angle)
    return Point(p.x*cos(angle)-p.y*sin(angle),p.x*sin(angle)+p.y*cos(angle))
end

function Base.:/(p1::Point, p2::Point)
    return Point(p1.x / p2.x, p1.y / p2.y)
end

function Base.:-(num::Number, p::Point)
    return Point(num - p.x, num - p.y)
end

function Base.:-(p1::Point, num::Number)
    return Point(p1.x - num, p1.y - num)
end

function Base.isless(num::Number, p::Point)
    return num < p.x && num < p.y
end

function Base.isless(p::Point, num::Number)
    return p.x < num && p.y < num
end

function Base.abs(p::Point)
    return Point(abs(p.x), abs(p.y))
end

function Base.isnan(p::Point)
    return isnan(p.x) || isnan(p.y)
end

function Base.min(p1::Point, p2::Point)
    return Point(min(p1.x, p2.x), min(p1.y, p2.y))
end

function Base.max(p1::Point, p2::Point)
    return Point(max(p1.x, p2.x), max(p1.y, p2.y))
end

function randF(min,max)
    return min+rand()*(max-min)
end

function randF(max)
    return randF(0,max)
end

function randI(min,max)
    return floor(randF(min,max))
end

function randI(max)
    return randI(0,max)
end

function linspace(start, stop, len)
    return range(start, stop=stop, length=len)
end

function meshgrid(x,y)
    return [Point(x,yi) for yi in y],[Point(xi,y) for xi in x]
end

function meshgrid(x_points::Int, y_points::Int, x_range::Tuple{Float64, Float64}, y_range::Tuple{Float64, Float64})
    x_lin = range(x_range[1], stop=x_range[2], length=x_points)
    y_lin = range(y_range[1], stop=y_range[2], length=y_points)
    
    x_grid = repeat(reshape(x_lin, 1, :), y_points, 1)
    y_grid = repeat(y_lin, 1, x_points)
    
    return x_grid, y_grid
end

function truncate(x)
    return x>0 ? floor(x) : ceil(x)
end

function Base.round(x,step)
    return round(x/step)*step
end

function parabola(x,root0=0,root1=1)
    a = 1 / ((0.5 * (root0 + root1) - root0) * (0.5 * (root0 + root1) - root1))
    return a*(x-root0)*(x-root1)/((root0-root1)^2)
end

function limit(x,x0,x1)
    return min(max(x,x0),x1)
end

function limitMag(x,x0)
    return limit(x,-x0,x0)
end

function norm(p::Point)
    return normalize(p)
end