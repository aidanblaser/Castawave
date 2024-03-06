using DrWatson
@quickactivate "Castawave"
using LaTeXStrings
using Printf
using MAT
using Trapz

include(projectdir()*"/src/MainSolver.jl")

n = 128
A = 0.3
Δt = 0.01
tf = 10.
L = 2π;
k = 1;
h = 0;
smoothing = false;

X = [(α * L / n) - A*sin(k*α*L/n) - A^3*k^2*sin(k*α*L/n) - A^4*k^3 / 3 * sin(2*k*α*L/n) for α in 1:n]
Y = [(cos(k * α*L / n )) * A + 1/6*A^4*k^3*cos(k*α*L/n) .+ 0*(A^2*k / 2).+ 0*A^4*k^3 * 1/2 for α in 1:n]
ϕ = [sqrt(GRAVITY/k) * A * exp.(k*Y[α]) * sin(k*X[α]) for α in 1:n]
# Use if you want a generic initial condition
# 
# using MAT 
if true
    vars = matread(projectdir()*"/data/ClamondICn128.mat")
    X = vec(vars["X"]);
    Y = vec(vars["Y"]);
    ϕ = vec(vars["F"]);
end

n = length(X)

scatter(X,Y)


MWL = sum(Y)/n
MWL = simpsons_rule_periodic(X,Y)/(2π)
sol= runSim(n, X, Y, ϕ, Δt, Float64(tf),L,h)



# END OF SIMULATION CODE (remaining code are tests or modified methods)

function visualize(interval::Int, fps::Int)
    anim = @animate for t = 0:Δt:tf
        xvals = sol(t)[1:n]
        yvals = sol(t)[n+1:2*n]
        ϕvals = sol(t)[2*n+1:3*n]

        scatter([xvals],[yvals], label = "Timestepped", legend = :bottomright,
         framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1,
          dpi = 300, xlabel=L"x \,(m)",ylabel=L"z \,(m)", title= @sprintf("Time: %.1f s", t),
          xlims=(0,L),ylims = (-1.9,1.9),aspect_ratio=1)
    end every interval
    gif(anim, projectdir()*"/plots/RK4noSmooth.gif", fps=fps)
end
visualize(1, 30)


energy, MWL_check = computeEnergy(sol,n,Δt,tf)
plot(0:Δt:tf,energy .- energy[1],xlabel=L"t \, (s)",ylabel =L"E - E_0 \quad (J/kg)",legend=false,title= L"E_0: \quad "*@sprintf("%.4f ", energy[1])*L"(J/kg)")
plot(0:Δt:tf,MWL_check,xlabel=L"t \, (s)",ylabel ="MWL"*L" \quad (m)",legend=false)
savefig(projectdir()*"/plots/Energy")