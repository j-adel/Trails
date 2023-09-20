
using Plots
include("fns.jl")

c = -5
x0=-2; x1=2; y0=3; y1=5
x_vals = linspace(x0,x1,100)
p = plot(x_vals, [scaledEase.(x_vals, c,x0,x1,y0,y1), scaledStep.(x_vals, c,x0,x1,y0,y1), 
                    modulo(x_vals, y0, y1), smoothStep.(x_vals)],
         label=["scaledEase" "scaledStep"],
         legend=:bottomright,
         lw=2)
# limit y range to [0,1]
title!("Easing Curves for c=2")
xlabel!("x")
ylabel!("Value")
plot!(p, xlims=(x0,x1),ylims=(y0,y1))
display(p)
