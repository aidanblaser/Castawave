using DrWatson
@quickactivate "Castawave"
using Plots
using JLD2
using LaTeXStrings

##### Retrieve wave data and kinematics
df = jldopen(projectdir()*"/data/a.jld2", "r")

df["sol"] = sol
df["n"] = n
df["A"] = A
df["Δt"] = Δt
df["tf"] = tf
df["L"] = L;
df["k"] = k;
df["h"] = h;
df["smoothing"] = smoothing;

x = sol(t)[1:n]
y = sol(t)[n+1:2*n]
ẋ = sol(t,Val{1})[1:n]
ẏ = sol(t,Val{1})[n+1:2*n]
ax = sol(t,Val{2})[1:n]
ay = sol(t,Val{2})[n+1:2*n]
ϕ = sol(t)[2*n+1:end]

##### Animate wave and highlight selected particles
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


##### Velocity and acceleration as function of time


##### Curve fitting



