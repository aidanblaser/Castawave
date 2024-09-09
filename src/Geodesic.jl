using DrWatson
@quickactivate "Castawave"
using Plots

include(projectdir()*"/src/MainSolver.jl")
include(projectdir()*"/src/ClamondIC.jl")
include(projectdir()*"/src/DoldMono.jl")

N = 256;
A = 0.4
X,Y,ϕ,c = getIC(Inf,A,N÷2)
tf = 1.6;
Δt = 1e-4
L = 2π
h = 1.2
tol = 1e-6
smoothed = true


scatter(X,Y)

Xfull, Yfull, ϕfull, t = @time runSim(N, X, Y, ϕ, Δt, Float64(tf),L,h,ϵ = tol,smoothing=smoothed)

gr()
scatter(Xfull[end,:],Yfull[end,:],background_color="black",aspect_ratio=1,
framestyle=:box,markerstrokewidth=0,markersize=1,dpi=300,xlabel="x (m)",ylabel="y(m)",legend=false)

plotlyjs()
plot(Xfull,t,Yfull .- 1/2*0 * GRAVITY *t,legend=false,zlims=(-5,5),xlabel="x (m)",ylabel= "t (s)",zlabel="y (m)")