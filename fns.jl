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

function randF(min,max)
    return min+rand()*(max-min)
end

function randF(max)
    return randF(0,max)
end
