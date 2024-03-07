using DrWatson
@quickactivate "Castawave"
using LaTeXStrings
using Printf
using JLD2

# This script runs the simulation for a wide range of parameter space, and saves the outputs
include(projectdir()*"/src/MainSolver.jl")

AmplitudeRange = [0.001,0.01,0.05,0.1,0.15,0.2]#,0.25,0.3]
L = 2π
n = 512
Δt = 0.01
tf = 10.
k = 1;
h = 0;

SolArray = [];
EArray = [];
MWL_Array = [];

for A ∈ AmplitudeRange
    println("Amplitude = $A")
    X = [(α * L / n) - A*sin(k*α*L/n) - A^3*k^2*sin(k*α*L/n) - A^4*k^3 / 3 * sin(2*k*α*L/n) for α in 1:n]
    Y = [(cos(k * α*L / n )) * A + 1/6*A^4*k^3*cos(k*α*L/n) .+ (A^2*k / 2).+ A^4*k^3 * 1/2 for α in 1:n]
    ϕ = [sqrt(GRAVITY/k) * A * exp.(k*Y[α]) * sin(k*X[α]) for α in 1:n]

    sol= runSim(n, X, Y, ϕ, Δt, Float64(tf),L,h)

    energy, MWL_check = computeEnergy(sol,n,Δt,tf)
    push!(SolArray,sol)
    push!(EArray,energy)
    push!(MWL_Array,MWL_check)
end

jldsave(projectdir()*"/data/CastawaveN512.jld2"; SolArray,EArray,MWL_Array)
