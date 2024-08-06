using DrWatson
@quickactivate "Castawave"
using LaTeXStrings
using Printf
using MAT
using Trapz
using Plots
plotlyjs()

using HDF5

include(projectdir()*"/src/MainSolver.jl")
include(projectdir()*"/src/ClamondIC.jl")
include(projectdir()*"/src/DoldMono.jl")

N = 1024
n = N
A = 0.3
k = 1;
tf = 5.0;
Δt = 1e-3
L = 2π
h = 0.0
tol = 1e-6
smoothed = true;

# X = [(α * L / n) - 0*A*sin(k*α*L/n) - 0*A^3*k^2*sin(k*α*L/n) - 0*A^4*k^3 / 3 * sin(2*k*α*L/n) for α in n÷2:n+n÷2-1]
# Y = [(cos(k * α*L / n )) * A + 0*1/6*A^4*k^3*cos(k*α*L/n) .+ 0*(A^2*k / 2).+ 0*A^4*k^3 * 1/2 for α in n÷2:n+n÷2-1]
# ϕ = [sqrt(9.81/k) * A   * sin(k*X[α]) for α in 1:n]
X = [(α * L / n) for α in n÷2:n+n÷2-1]
Y = [(cos(k * α*L / n )) * A for α in n÷2:n+n÷2-1]
ϕ = [sqrt(9.81/k) * A * sin(k*X[α]) for α in 1:n]


Xfull, Yfull, ϕfull, t = runSim(N, X, Y, ϕ, Δt, Float64(tf),L,h,ϵ = tol,smoothing=smoothed)
scatter(Xfull[end,:], Yfull[end,:])

for i in 1:101
    A = 0.3 + ((i-1) / 1000)
    X = [(α * L / n) for α in n÷2:n+n÷2-1]
    Y = [(cos(k * α*L / n )) * A for α in n÷2:n+n÷2-1]
    ϕ = [sqrt(9.81/k) * A * sin(k*X[α]) for α in 1:n]

    Xfull, Yfull, ϕfull, t = runSim(N, X, Y, ϕ, Δt, Float64(tf),L,h,ϵ = tol,smoothing=smoothed)

    fid = h5open("lw0_"*string(Int64(299+i))*".h5","w")
    fid["A"] = A
    fid["N"] = N
    fid["X"] = Xfull
    fid["Y"] = Yfull
    fid["Phi"] = ϕfull
    fid["t"] = t
    close(fid)
end
