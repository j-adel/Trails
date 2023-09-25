using Plots

# Define the ODE system
dxdt(x, y) = y
dydt(x, y) = -x

function leapfrog_integration(x0, y0, T, dt)
    x = x0
    y = y0
    xs, ys = [x0], [y0]
    times = 0:dt:T

    # # Initial half step for velocity (Here, y acts as a velocity for x and vice versa)
    # y += 0.5 * dt * dydt(x, y)
    # x += 0.5 * dt * dxdt(x, y)

    for t in times[2:end]
        # Update positions
        x += dt * dxdt(x, y)
        y += dt * dydt(x, y)

        
        # Update velocities using the new positions
        y += dt * dydt(x, y)
        x += dt * dxdt(x, y)
        push!(xs, x)
        push!(ys, y)
    end

    return xs, ys
end

function eulerIntegration(x0, y0, T, dt)
    x = x0
    y = y0
    xs, ys = [x0], [y0]
    times = 0:dt:T

    for (i, t) in enumerate(times[2:end])
        # Update positions
        x += dt * dxdt(xs[end], ys[end])
        y += dt * dydt(xs[end], ys[end])
        println("r ",sqrt(x^2+y^2))
        push!(xs, x)
        push!(ys, y)
    end

    return xs, ys
end

# Initial conditions
x0, y0 = .71, .71

# Time span and time step
T = 10.0
dt = .1

xs, ys = leapfrog_integration(x0, y0, T, dt)
xe, ye = eulerIntegration(x0, y0, T, dt)
plot(xs, ys, label="Leapfrog trajectory", xlabel="x", ylabel="y", legend=:topright)
plot!(xe, ye, label="Euler trajectory", xlabel="x", ylabel="y", legend=:topright)
# make aspec ratio 1:1
plot!(aspect_ratio=:equal)