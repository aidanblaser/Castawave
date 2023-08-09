using DrWatson
@quickactivate "PacketDrift"
using LaTeXStrings
using Printf


include(projectdir()*"/src/MainSolver.jl")

n = 256
A = 0.01
Δt = 0.01
tf = 2
L = 2π;
k = 1;

X = [(α * 2π / n) - A*sin(k*α*2π/n) - A^3*k^2*sin(k*α*2π/n) - A^4*k^3 / 3 * sin(2*k*α*2π/n) for α in 1:(2*n)]
Y = [(cos(k * α*2π / n )) * A + 1/6*A^4*k^3*cos(k*α*2π/n) + A^2*k / 2 + A^4*k^3 / 2 for α in 1:(2*n)]
ϕ = [sqrt(GRAVITY/k) * A * exp.(k*Y[α]) * sin(k*X[α]) for α in 1:(2*n)]

xf, yf, ϕf, wl, ta = run(n*2, X, Y, ϕ, Δt, Float64(tf),2*L)
jldsave("RK4.1.jld2"; x=xf, y=yf, ϕ=ϕf, N=n, A=A, dt=Δt, tf=tf)


function visualize(interval::Int, fps::Int,tf,dt)
    anim = @animate for i ∈ 1:(ceil(Int, tf/dt) + 2)
        # scatter([sol[i][:,1]], [sol[i][:,2]], label = "Timestepped", legend = :bottomright, framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1, dpi = 300, xlabel=L"x \,(m)",ylabel=L"\eta \,(m)", title= @sprintf("Time: %.3f s", (i-1)*dt))

        scatter([xf[i,:]], [yf[i,:]], label = "Timestepped", legend = :bottomright, framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1, dpi = 300, xlabel=L"x \,(m)",ylabel=L"z \,(m)", title= @sprintf("Time: %.3f s", (i-1)*Δt))
        scatter!([xf[1,:]], [yf[1,:]], label = "Initial position", framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1, dpi = 300, xlabel=L"x \,(m)",ylabel=L"\eta \,(m)", title= @sprintf("Time: %.3f s", (i-1)*Δt))
    end every interval
    gif(anim, "RK4.1.gif", fps=fps)
end
visualize(10, 5,tf,Δt)

# END OF SIMULATION CODE (remaining code are tests or modified methods)