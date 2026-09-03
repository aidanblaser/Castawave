#=
QuickStart.jl

The fastest path from a clean checkout to a running simulation: load the
Castawave package, build a small analytic initial wave, run it, and take
a quick look at the result. Start here if you're new to this codebase, or
just want to sanity-check that your Julia environment is set up correctly.

This script lives in its own nested environment (scripts/Project.toml),
which already has Castawave set up as a dev dependency plus Plots for the
visualization at the end. The first time you use it, from the repo root:

    julia --project=scripts -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

After that, run it either as a script:
    julia --project=scripts scripts/QuickStart.jl
or paste it into a REPL:
    julia --project=scripts
    julia> include("scripts/QuickStart.jl")
=#

using Castawave
using Plots

# ---------------------------------------------------------------------
# 1. Set up the run parameters
# ---------------------------------------------------------------------
N  = 256    # number of surface particles (higher N = finer spatial resolution)
L  = 2π     # domain length / spatial periodicity, in meters
h  = 0.0    # still-water depth; 0.0 means infinite depth
A  = 0.1    # initial wave amplitude, in meters
k  = 1      # initial wavenumber
Δt = 0.01   # interval, in seconds, at which output is saved. NOT the
            # internal integration step size - that's set adaptively from
            # errortol below. See SimulationParameters' docstring in
            # src/MainSolver.jl for the full explanation.
T  = 5.0    # duration to simulate, in seconds

p = SimulationParameters(L, h, Δt, T, smoothing=true)

# ---------------------------------------------------------------------
# 2. Build an initial condition
# ---------------------------------------------------------------------
# A small-amplitude, approximately-Stokes traveling wave (the same
# analytic form used in scripts/BasicSolver.jl). Particle labels ξ are
# just the evenly-spaced indices 1:N; X, Y and ϕ are each length-N vectors
# giving that particle's initial horizontal position, vertical position,
# and velocity potential.
X = [(α * L / N) - A*sin(k*α*L/N) for α in 1:N]
Y = [A*cos(k*α*L/N) for α in 1:N]
ϕ = [sqrt(9.81/k) * A * sin(k*α*L/N) for α in 1:N]

# ---------------------------------------------------------------------
# 3. Run the simulation
# ---------------------------------------------------------------------
Xfull, Yfull, ϕfull, t = runSim(X, Y, ϕ, p)

println("Ran to t = ", round(t[end], digits=3), " s, saving ", length(t), " timesteps.")

# ---------------------------------------------------------------------
# 4. Take a quick look
# ---------------------------------------------------------------------
plt = plot(Xfull[1, :], Yfull[1, :], label="t = 0",
           xlabel="x (m)", ylabel="y (m)", aspect_ratio=:equal,
           title="Castawave quick start")
plot!(plt, Xfull[end, :], Yfull[end, :], label="t = $(round(t[end], digits=2)) s")
display(plt)
