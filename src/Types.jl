#=
This file defines SimulationParameters (the user-facing run configuration
struct) and SolverWorkspace (a bundle of preallocated scratch arrays and
per-run constants that get threaded through the timestepping hot path:
fixedTimeOperations!, TimeDerivatives!, ABMatrices!, NormalInversion!,
RealPhi!). Previously each of those functions took its own long list of
loose positional arguments for these buffers (30-40 in the worst cases);
bundling them here means each function just takes `ws` plus the few
arguments that actually vary between calls (the current X/Y/ϕ vs. the
predicted X/Y/ϕ, and which output buffers to write into).

All fields are concretely typed, so accessing ws.foo compiles to the same
machine code as a plain local variable - this struct is purely an
organizational convenience over the preallocation pattern already in use,
not a new source of allocation or indirection. Included before
HelperFunctions.jl so that every function that annotates an argument as
::SolverWorkspace already has the type available.
=#

# Create a structure for the parameters so it knows it's invariant
struct SimulationParameters

    # Physical parameters
    L::Float64 # Length of domain in meters (spatial periodicity)
    h::Float64 # Depth of fluid at rest in meters
    dt::Float64 # Interval of physical time (seconds) between saved output snapshots. Does NOT bound the internal step size; see atol below and runSim's docstring.
    T::Float64 # Desired duration of simulation in seconds (may abort early if wave breaks)

    # Parameters re-scaled to have L=2π
    lengthScale::Float64
    timeScale::Float64
    h̃::Float64 
    dt̃::Float64 
    T̃::Float64

    # Parameters with default values 
    atol::Float64 # absolute tolerance on predictor/corrector disagreement in X/Y/ϕ - sets the internal timestep; see runSim's docstring
    errortol::Float64 # no longer sets the internal timestep - only gates conditional smoothing (via erm)
    smoothing::Bool 
    g::Float64
    stabilityFactor::Float64 # Safety factor, roughly 0 to 1, on the Dold (1992, JCP 103, 90-115)
                              # strong-instability timestep threshold. 0.6 matches
                              # Dold's own recommended (most conservative) setting;
                              # raise it (up to ~1.2) only if you've confirmed runs
                              # stay stable and want to take larger steps.
end

function SimulationParameters(L,h,dt,T;atol=1e-6,errortol=1e-4,smoothing=true,g=9.81,stabilityFactor=0.6)
        # Convert variables to explicit types 
        L = Float64(L)
        h = Float64(h)
        dt = Float64(dt)
        T = Float64(T)
        atol = Float64(atol)
        errortol = Float64(errortol)
        smoothing = Bool(smoothing)
        g = Float64(g)
        stabilityFactor = Float64(stabilityFactor)

        # Parameters scaled to have 2π periodicity
        lengthScale = Float64(L/2π)
        timeScale = sqrt(lengthScale)
        h̃ = h/lengthScale
        dt̃ = dt/timeScale
        T̃ = T/timeScale
        return SimulationParameters(L, h, dt, T,
                               lengthScale, timeScale, h̃, dt̃, T̃,
                               atol, errortol, smoothing, g, stabilityFactor)
end

struct SolverWorkspace
    # --- Fixed configuration for the run (never reassigned after construction) ---
    N::Int
    H::Float64
    g::Float64
    c1::Vector{Float64}              # 1st-derivative (5-pt) stencil coefficients
    c2::Vector{Float64}              # 2nd-derivative (6-pt) stencil coefficients
    smoothCoefficients::Vector{Float64}   # profile/potential smoothing stencil (Dold's m=15 formula)
    roughnessCoefficients::Vector{Float64} # fixed 11-point (m=11) formula used only for
                                            # the roughness diagnostic, independent of
                                            # whatever stencil smoothCoefficients uses
    spreadCoefficients::Vector{Float64}    # all-positive 11-point (m=11) windowing formula
                                            # (Dold's `spreadd`) used to spatially spread
                                            # the per-point roughness weight for smthwd-style
                                            # weighted smoothing - distinct from the signed
                                            # roughnessCoefficients stencil above

    # --- Conformally mapped geometry (complex) ---
    Ω::Vector{ComplexF64}
    Ω_ξ::Vector{ComplexF64}
    Ω_ξξ::Vector{ComplexF64}

    # --- Real-plane geometry derivatives ---
    X_ξ::Vector{Float64}
    Y_ξ::Vector{Float64}
    inverseds2::Vector{Float64}

    # --- Dense BEM matrices, rebuilt every call by ABMatrices! ---
    A::Matrix{Float64}
    B::Matrix{Float64}
    Cinter::Matrix{Float64}
    C::Matrix{ComplexF64}
    ΔΩ::Matrix{ComplexF64}

    # --- Normal-derivative inversion scratch (NormalInversion!) ---
    ϕ_ξ::Vector{Float64}
    ϕ_ξξ::Vector{Float64}
    ϕ_ν::Vector{Float64}
    b::Vector{Float64}

    # --- TimeDerivatives! scratch: Eulerian time derivatives & velocity gradients ---
    ϕ_t::Vector{Float64}
    ϕ_tξ::Vector{Float64}
    ϕ_tξξ::Vector{Float64}
    ϕ_tν::Vector{Float64}
    u_t::Vector{Float64}
    u_tξ::Vector{Float64}
    v_t::Vector{Float64}
    v_tξ::Vector{Float64}
    u_ξ::Vector{Float64}
    v_ξ::Vector{Float64}
    u_x::Vector{Float64}
    v_x::Vector{Float64}
    u_tx::Vector{Float64}
    v_tx::Vector{Float64}
    u_xξ::Vector{Float64}
    v_xξ::Vector{Float64}
    u_xx::Vector{Float64}
    v_xx::Vector{Float64}
    ϕ_tt::Vector{Float64}
    ϕ_ttξ::Vector{Float64}
    ϕ_ttξξ::Vector{Float64}
    ϕ_ttν::Vector{Float64}
    u_tt::Vector{Float64}
    v_tt::Vector{Float64}

    # --- 11-point smoothing scratch ---
    Ω_sm_temp::Vector{Float64}

    # --- smthwd (weighted smoothing) scratch ---
    roughnessVec::Vector{Float64}          # per-point roughness r(i) from roughnessVector!
    roughnessWeightTemp::Vector{Float64}   # scratch: clipped/normalized roughness before spreadd
    roughnessWeight::Vector{Float64}       # final per-point smoothing weight h(i) in [0,1]
end

function SolverWorkspace(N::Int, H::Real, g::Real)
    r() = Vector{Float64}(undef, N)
    z() = Vector{ComplexF64}(undef, N)

    c1 = [2100.0, -600.0, 150.0, -25.0, 2.0] ./ 2520.0
    c2 = [-9220.75, 5250.0, -750.0, 125.0, -15.625, 1.0] ./ 3150.0
    smoothCoefficients = [3432.0, -3003.0, 2002.0, -1001.0, 364.0, -91.0, 14.0, -1.0] ./ (2^14)
    # Dold's `rough` diagnostic always uses an 11-point formula (m=11),
    # regardless of what m the profile/potential smoothing above uses.
    roughnessCoefficients = [252.0, -210.0, 120.0, -45.0, 10.0, -1.0] ./ 1024.0
    # Same magnitude coefficients as roughnessCoefficients, but all-positive:
    # this is a plain windowed average (Dold's `spreadd`), not a signed residual.
    spreadCoefficients = [252.0, 210.0, 120.0, 45.0, 10.0, 1.0] ./ 1024.0

    return SolverWorkspace(
        N, Float64(H), Float64(g), c1, c2, smoothCoefficients, roughnessCoefficients, spreadCoefficients,
        # Ω, Ω_ξ, Ω_ξξ
        z(), z(), z(),
        # X_ξ, Y_ξ, inverseds2
        r(), r(), r(),
        # A, B, Cinter, C, ΔΩ
        Matrix{Float64}(undef, N, N), Matrix{Float64}(undef, N, N), Matrix{Float64}(undef, N, N),
        Matrix{ComplexF64}(undef, N, N), Matrix{ComplexF64}(undef, N, N),
        # ϕ_ξ, ϕ_ξξ, ϕ_ν, b
        r(), r(), r(), r(),
        # ϕ_t, ϕ_tξ, ϕ_tξξ, ϕ_tν
        r(), r(), r(), r(),
        # u_t, u_tξ, v_t, v_tξ
        r(), r(), r(), r(),
        # u_ξ, v_ξ, u_x, v_x
        r(), r(), r(), r(),
        # u_tx, v_tx
        r(), r(),
        # u_xξ, v_xξ
        r(), r(),
        # u_xx, v_xx
        r(), r(),
        # ϕ_tt, ϕ_ttξ, ϕ_ttξξ, ϕ_ttν
        r(), r(), r(), r(),
        # u_tt, v_tt
        r(), r(),
        # Ω_sm_temp
        r(),
        # roughnessVec, roughnessWeightTemp, roughnessWeight
        r(), r(), r(),
    )
end
