using DrWatson
@quickactivate "Castawave"
using LaTeXStrings
using Printf


include(projectdir()*"/src/MainSolver.jl")

n = 512
A = 0.5
Δt = 0.01
tf = 2.84
L = 2π;
k = 1;
h = 0;  #infinite depth

X = [(α * L / n) - A*sin(k*α*L/n) - A^3*k^2*sin(k*α*L/n) - A^4*k^3 / 3 * sin(2*k*α*L/n) for α in 1:n]
Y = [(cos(k * α*L / n )) * A + 1/6*A^4*k^3*cos(k*α*L/n) + A^2*k / 2 + A^4*k^3 / 2 for α in 1:n]
ϕ = [sqrt(GRAVITY/k) * A * exp.(k*Y[α]) * sin(k*X[α]) for α in 1:n]

using MAT 
vars = matread(projectdir()*"/data/ClamondIC.mat")
X = vec(vars["X"]);
Y = vec(vars["Y"]);
ϕ = vec(vars["F"]);


xf, yf, ϕf,time, wl, ta = runsim(n, X, Y, ϕ, Δt, Float64(tf),L,h)
jldsave(projectdir()*"/data/RK4.1.jld2"; x=xf, y=yf, ϕ=ϕf, N=n, A=A, dt=Δt, tf=tf)


function visualize(interval::Int, fps::Int,time)
    anim = @animate for i ∈ 1:length(time)
        # scatter([sol[i][:,1]], [sol[i][:,2]], label = "Timestepped", legend = :bottomright, framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1, dpi = 300, xlabel=L"x \,(m)",ylabel=L"\eta \,(m)", title= @sprintf("Time: %.3f s", (i-1)*dt))

        scatter([xf[i,:]], [yf[i,:]], label = "Timestepped", legend = :bottomright, framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1, dpi = 300, xlabel=L"x \,(m)",ylabel=L"z \,(m)", title= @sprintf("Time: %.3f s", time[i]))
        scatter!([xf[1,:]], [yf[1,:]], label = "Initial position", framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1, dpi = 300, xlabel=L"x \,(m)",ylabel=L"\eta \,(m)", title= @sprintf("Time: %.3f s", time[i]))
    end every interval
    gif(anim, projectdir()*"/plots/RK4.1.gif", fps=fps)
end
visualize(10, 5,time)

# END OF SIMULATION CODE (remaining code are tests or modified methods)