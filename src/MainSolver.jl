#=
This file contains the main routine for initialising and running the adapted Dold flow solver. The program is initialized through the input parameters of the run function, which then runs the solver for the specified time or amount.

Package dependencies as well as the other two files needed to run the solver are initially called and included. The fixedTimeOperations computes all necessary quantities needed to step the model forward in time, while the overarching run function marches the X, Y, and ϕ Vectors in time until outputting a matrix wherein each row corresponds to a fixed time. 

Note that, during initialization, particle labels are taken as the indices of the X, Y, and ϕ vectors, and are as such evenly spaced no matter the positional values.
=#

using DrWatson
@quickactivate "Castawave"
using LinearAlgebra

include(projectdir()*"/src/Constants.jl")
include(projectdir()*"/src/Types.jl")
include(projectdir()*"/src/HelperFunctions.jl")

#= TODO 
- Better implement time stepping so it lands on intervals of dt (even if it takes steps between)
- 
=#

# SimulationParameters now lives in Types.jl, alongside SolverWorkspace.
# It has to be defined before HelperFunctions.jl is included below, because
# a couple of legacy (non-mutating) functions there - PhiTimeDer and
# TimeDerivatives - annotate an argument as ::SimulationParameters, and
# Julia resolves type annotations in a method signature at the time that
# `function` definition is evaluated, not lazily. Defining the struct here
# instead (after HelperFunctions.jl is already included) throws
# `UndefVarError: SimulationParameters not defined` from inside
# HelperFunctions.jl on a fresh run.


function fixedTimeOperations(X::AbstractVector{<:Real}, Y::AbstractVector{<:Real}, ϕ::AbstractVector{<:Real}, p::SimulationParameters,N::Int,H)
    #=
    The fixedTimeOperations function wraps all of the necessary operations for finding the quantities needed for the next timestep. These consist of the finding the R_ξ derivative and the ϕ_ξ, ϕ_ν derivatives, from which ϕ_t is given from Bernoulli's condition. Then, the change in both R = X + iY and ϕ is known and the system can be evolved to the next timestep.
    
    It has many calls to the underlying helper functions file for clarity and compartementalization of the code.
    
    Input:
    X - real vector of initial particle x-positions on the surface
    Y - real vector of initial particle y-positions on the surface
    ϕ - real vector of scalar velocity potential for particles on the surface
    p - SimulationParameters structure

    Output:
    ϕ_x - real vector of U velocity for particles on the surface
    ϕ_y - real vector of V velocity for particles on the surface
    ϕ_D - real vector of the material derivative of ϕ from the dynamic boundary condition
    (higher derivatives of X,Y,ϕ are outputted as well)
    =#

    # Compute Ω by taking a conformal map 
    conformalMap!(Ω,X .+ im*Y)

    # Compute derivatives 
    # Requires to be non-dimensionsionalized
    R_ξ = DDI1(X .+ im*Y,2π)
    Ω_ξ = DDI1(Ω,0.0)
    Ω_ξξ = DDI2(Ω,0.0)


    # The matrix method described in Dold is used to find the normal derivative of the potential.
    A, B, ℵ = ABMatrices(Ω, Ω_ξ, Ω_ξξ, H)
    ϕ_ξ, ϕ_ν = NormalInversion(ϕ, A, ℵ)

    # All necessary information having been computed, we transform back to the real frame to output conveniant timestepping quantities.
    ϕ_x, ϕ_y = RealPhi(R_ξ, ϕ_ξ, ϕ_ν)
    
    # Use all this to compute up to third order time derivatives
    TimeDerivatives!(DϕDt,DuDt,DvDt,D2ϕDt2,D2uDt2,D2vDt2,D3ϕDt3,
    X_ξ,Y_ξ,ϕ_x,ϕ_y,A,b,ℵ,Y,p,N,c1,c2,
    ϕ_t,ϕ_tξ,ϕ_tν,u_t,v_t,u_ξ,v_ξ,u_x,v_x,u_tx,v_tx,u_xξ,v_xξ,u_xx,v_xx,
    ϕ_tt,ϕ_ttξ,ϕ_ttν,u_tt,v_tt)
end

function fixedTimeOperations!(ws::SolverWorkspace, X, Y, ϕ, ϕ_x, ϕ_y,
    DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3)
    #=
    The fixedTimeOperations function wraps all of the necessary operations for finding the quantities needed for the next timestep. These consist of the finding the R_ξ derivative and the ϕ_ξ, ϕ_ν derivatives, from which ϕ_t is given from Bernoulli's condition. Then, the change in both R = X + iY and ϕ is known and the system can be evolved to the next timestep.
    
    It has many calls to the underlying helper functions file for clarity and compartementalization of the code.
    
    Input:
    ws - SolverWorkspace holding every preallocated scratch buffer plus N, H, g, c1, c2 for this run
    X - real vector of particle x-positions on the surface (current or predicted)
    Y - real vector of particle y-positions on the surface (current or predicted)
    ϕ - real vector of scalar velocity potential for particles on the surface (current or predicted)

    Output (written in place into the buffers passed in):
    ϕ_x - real vector of U velocity for particles on the surface
    ϕ_y - real vector of V velocity for particles on the surface
    DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3 - material derivatives up to third order
    =#
    (; Ω, Ω_ξ, Ω_ξξ, X_ξ, Y_ξ, inverseds2, N, c1, c2, ϕ_ξ, ϕ_ξξ, ϕ_ν) = ws

    # Compute Ω by taking a conformal map 
    conformalMap!(Ω,X .+ im*Y)
    # Compute derivatives 
    DDI1!(X_ξ,X,2π,N,c1)
    DDI1!(Y_ξ,Y,0.0,N,c1)
    DDI1!(Ω_ξ,Ω,zero(eltype(Ω)),N,c1)
    DDI2!(Ω_ξξ,Ω,zero(eltype(Ω)),N,c2)

    # Computing inverse spacing (used a lot)
    @inbounds for i ∈ 1:N 
        inverseds2[i] = (X_ξ[i]^2 + Y_ξ[i]^2)^(-1)
    end 

    # The matrix method described in Dold is used to find the normal derivative of the potential.
    ℵ = ABMatrices!(ws, Ω, Ω_ξ, Ω_ξξ)
    NormalInversion!(ws, ϕ_ξ, ϕ_ξξ, ϕ_ν, ϕ, ℵ)

    # All necessary information having been computed, we transform back to the real frame to output conveniant timestepping quantities.
    RealPhi!(ws, ϕ_x, ϕ_y, ϕ_ξ, ϕ_ν)

    # Compute up to third order time derivatives 
    TimeDerivatives!(ws, ℵ, ϕ_x, ϕ_y, Y, DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3)
end


function runSim(X::AbstractVector{<:Real}, Y::AbstractVector{<:Real}, ϕ::AbstractVector{<:Real}, p::SimulationParameters)
    #=
    The run function is the master function of the program, taking in the initial conditions for the system
    and timestepping it forward until some final time. The outputs are written into arrays
    (which can also be exported through JDL2 for larger files).
    
    Input:
    X  - real vector (length N) of initial particle x-positions on the surface
    Y  - real vector (length N) of initial particle y-positions on the surface
    ϕ  - real vector (length N) of initial scalar velocity potential for particles on the surface
    p  - Simulation parameters structure. Includes the fields:
        L - periodicity of domain (meters)
        h - depth of fluid at rest (meters)
        dt - interval of physical time (seconds) at which the output series is
             saved. Does NOT set the internal integration step size any more -
             that's governed entirely by errortol (and the stability cap
             below); dt only controls how densely the returned series is
             sampled. Internal steps are taken freely between save points and
             simply aren't recorded, with the final sub-step of each interval
             shortened so the saved series lands on clean multiples of dt.
        tf - duration of simulation (seconds, will abort early if wave breaks)
        errortol - error tolerance, defaults to 1e-5. Lower values mean greater accuracy. Sets the actual internal timestep based on nonlinearity
        smoothing - Default to true. At each timestep, applies an 11-pt smoothing filter to remove "sawtooth" modes which tend to appear in these codes 
        g - value of gravitational acceleration (defaults to 9.81)

    Output:
    X_timeseries - real array of x-positions for particles on the surface at each time step
    X_timeseries - real array of y-positions for particles on the surface at each time step
    ϕ_timeseries - real vector of ϕ values for particles on the surface at each time step
    time         - real array of timesteps from 0 to tf with step size dt
    =#
    N = length(X);

    # For simplicity, non-dimensionalize the data 
    lengthScale = p.lengthScale
    timeScale = p.timeScale
    XS = X / lengthScale
    YS = Y / lengthScale 
    ϕS = ϕ / lengthScale / timeScale
    hS = p.h / lengthScale 
    H = conformalDepth(hS)
    gravity = p.g

    # Add breaking parameter 
    breaking = false

    # Initialize time vector
    t = [0.0]

    # Shift Y so that mean water level is 0 
    MWL = sum(YS .* DDI1(XS,2π))/N 
    
    # Create and initialize the timeseries fields (the *saved* output
    # series - sparse, spaced dt̃val apart in physical time; see Xcur/Ycur/
    # ϕcur below for the state that actually gets advanced every internal
    # step). 
    Xfull = Vector{Vector{Float64}}()
    Yfull = Vector{Vector{Float64}}()
    ϕfull = Vector{Vector{Float64}}()
    push!(Xfull,XS)
    push!(Yfull,YS .- MWL)
    push!(ϕfull,ϕS)

    # The actual "current" surface state, advanced every accepted internal
    # step regardless of whether that step lands on a save point. Kept
    # separate (and copied, not aliased) from Xfull/Yfull/ϕfull so that
    # internal sub-steps between save points don't get recorded.
    Xcur = copy(XS); Ycur = copy(YS .- MWL); ϕcur = copy(ϕS)
    tcur = 0.0
    # nextSaveTime is set from dt̃val a bit further down, once p.dt̃ has
    # actually been read into dt̃val.
    # prevΔt tracks the previous accepted step size, for the growth-rate
    # cap below (Dold's nicedt: dt = 0.d0 initially, so his growth cap is
    # likewise skipped on the very first step).
    prevΔt = 0.0

    # Preallocate every scratch buffer used by the timestepper into one workspace
    # (previously ~65 separate `similar(X)`/`Vector{...}(undef,N)` locals threaded
    # through fixedTimeOperations!/TimeDerivatives!/etc. as loose positional args).
    ws = SolverWorkspace(N, H, gravity)

    # These stay outside the workspace because both the "current" and "predicted"
    # sets need to exist simultaneously across the predictor-corrector step.
    Xnext = similar(X); Ynext = similar(X); ϕnext = similar(X)
    Xcorr = similar(X); Ycorr = similar(X); ϕcorr = similar(X)

    ϕ_x = similar(X); ϕ_y = similar(X)
    DϕDt = similar(X); DuDt = similar(X); DvDt = similar(X)
    D2ϕDt2 = similar(X); D2uDt2 = similar(X); D2vDt2 = similar(X); D3ϕDt3 = similar(X)

    ϕ_xp = similar(X); ϕ_yp = similar(X)
    DϕDtp = similar(X); DuDtp = similar(X); DvDtp = similar(X)
    D2ϕDt2p = similar(X); D2uDt2p = similar(X); D2vDt2p = similar(X); D3ϕDt3p = similar(X)

    T̃val = p.T̃
    smoothingval = p.smoothing 
    errortolval = p.errortol
    dt̃val = p.dt̃
    stabilityFactorVal = p.stabilityFactor
    nextSaveTime = dt̃val

    # ABMatrices! (called every timestep below) now parallelizes its O(N^2)
    # matrix build with Threads.@threads, and the lu! factorization that
    # follows it is BLAS/LAPACK-backed and can itself use multiple threads.
    # Both require Julia to actually have been started with more than one
    # thread (`julia -t auto`, or the JULIA_NUM_THREADS env var) - this is a
    # one-time check to surface it if that wasn't done, since it's an easy
    # free speedup to miss silently.
    if Threads.nthreads() == 1
        @info "Julia was started with only 1 thread (Threads.nthreads() == 1). ABMatrices! can use multiple threads for its per-timestep matrix build - restart with `julia -t auto` (or set JULIA_NUM_THREADS) to take advantage of it."
    elseif LinearAlgebra.BLAS.get_num_threads() == 1
        @info "LinearAlgebra.BLAS.get_num_threads() == 1: the per-timestep lu! factorization in ABMatrices! is BLAS/LAPACK-backed and can use multiple threads automatically. Consider calling `LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())` before running for a likely free speedup on that solve."
    end

    while tcur <= T̃val && !breaking
        try
            # Compute up to third order derivatives of X, Y, ϕ 
            fixedTimeOperations!(ws, Xcur, Ycur, ϕcur, ϕ_x, ϕ_y,
            DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3)

            # Roughness diagnostic (Dold's "rough"/erm), used below to decide
            # whether the surface needs smoothing this step. Always computed
            # from the fixed-order-3 accelerations, using the dedicated
            # 11-point (m=11) roughness stencil (independent of whatever
            # order the actual smoothing filter below uses).
            erm = maxRoughness(D2uDt2, D2vDt2, D3ϕDt3, N, ws.roughnessCoefficients)

            # Determing Adaptive Timestep. dt̃val is no longer a cap here - it's
            # the save interval, applied below via the landing clip instead
            # (which bounds Δt by at most dt̃val anyway, since it never lets
            # a step overshoot the next save point).
            thirdOrderMax = max(maximum(abs, D2uDt2),maximum(abs, D2vDt2),maximum(abs, D3ϕDt3))
            Δt = (errortolval*factorial(3)/thirdOrderMax)^(1/3)

            # Numerical-stability threshold on the timestep, ported from Dold's
            # "timstp"/strong-instability check (Dold 1992, J. Comp. Phys. 103,
            # 90-115). Explicit surface-tracking schemes like this one become
            # strongly unstable above a step size set by the local balance of
            # surface curvature/acceleration against gravity - independent of
            # the accuracy-based Δt above, which can look fine while a run is
            # quietly heading into instability (especially as a wave steepens
            # toward breaking). stabilityFactorVal defaults to 0.6, matching
            # Dold's own recommended (most conservative) setting.
            stabilityMax = 0.0
            @inbounds for i ∈ 1:N
                tp = abs(ws.X_ξ[i]*(DvDt[i]+gravity) - ws.Y_ξ[i]*DuDt[i]) / (ws.X_ξ[i]^2 + ws.Y_ξ[i]^2)
                stabilityMax = max(stabilityMax, tp)
            end
            Δt = min(Δt, stabilityFactorVal / sqrt(stabilityMax))

            # Growth-rate cap, ported from Dold's nicedt: "restricting the
            # rate of time-step growth to quadruple in about 5 steps" - Δt
            # is never allowed to grow by more than 32% over the previous
            # accepted step, regardless of what the error-tolerance/
            # stability formulas above would otherwise now allow. This is
            # what keeps the step size from snapping straight back up right
            # after a rough patch (e.g. near breaking) even once the
            # instantaneous curvature/acceleration briefly relaxes - it has
            # to ramp back up gradually instead. Skipped on the very first
            # step (prevΔt starts at 0, matching Dold's dt = 0.d0 init).
            if prevΔt > 0.0
                Δt = min(Δt, 1.32*prevΔt)
            end

            # Minimum timestep 1e-4
            Δt = max(Δt,1e-4)

            # Clip the step so it approaches the next save-time (dt̃val-
            # spaced) boundary. This matches Dold's own "pts" printout
            # interval explicitly - a fixed interval (|pts| > 0 "demands
            # some time-steps to land on time multiples of |pts|"), with no
            # breaking-triggered densification (not a feature of his code;
            # his only variant is a negative-pts mode adding extra prints
            # at kinetic-energy extrema, unrelated to breaking) - including
            # his nicedt "homing in smoothly from as many as almost 4 steps
            # distant": rather than taking full-size steps right up to the
            # boundary and clipping one abrupt final step, the step size
            # eases down over the last few steps (quartered, thirded,
            # halved, then exact) as it approaches. Internal steps between
            # save points are taken freely (per errortolval and the
            # stability/growth caps above) and simply not recorded. Only
            # home in on a save point that's still within the requested
            # duration; past the last full interval before T̃val, steps
            # proceed unclipped until tcur exceeds T̃val and the outer loop
            # stops.
            landingStep = false
            if nextSaveTime <= T̃val
                if tcur + 1.1*Δt >= nextSaveTime
                    Δt = nextSaveTime - tcur
                    landingStep = true
                elseif tcur + 1.9*Δt >= nextSaveTime
                    Δt = 0.5*(nextSaveTime - tcur)
                elseif tcur + 2.9*Δt >= nextSaveTime
                    Δt = (1.0/3.0)*(nextSaveTime - tcur)
                elseif tcur + 3.9*Δt >= nextSaveTime
                    Δt = 0.25*(nextSaveTime - tcur)
                end
            end

            # Use derivatives up to third order to make a predictor step 
            @inbounds for i ∈ 1:N
                Xnext[i] = Xcur[i] .+ ϕ_x[i] * Δt .+ Δt^2/2*DuDt[i] .+ Δt^3/6*D2uDt2[i] 
                Ynext[i] = Ycur[i] .+ Δt * (ϕ_y[i]) .+ Δt^2/2*DvDt[i] .+Δt^3/6*D2vDt2[i]
                ϕnext[i] = ϕcur[i] .+ Δt * (DϕDt[i]) .+ Δt^2/2*D2ϕDt2[i] .+Δt^3/6*D3ϕDt3[i]
            end

            # Estimate derivatives at the predicted surface
            fixedTimeOperations!(ws, Xnext, Ynext, ϕnext, ϕ_xp, ϕ_yp,
            DϕDtp, DuDtp, DvDtp, D2ϕDt2p, D2uDt2p, D2vDt2p, D3ϕDt3p)

            # Use predictor-corrector to average derivatives at predicted surface and current surface (like trapezoidal rule)
            @inbounds for i ∈ 1:N
                Xcorr[i] = Xcur[i] .+ Δt/2 *(ϕ_x[i] .+ ϕ_xp[i]) .+ Δt^2 / 12 *(DuDt[i] .- DuDtp[i]) .+ Δt^3 /24 * (D2uDt2[i] .+ D2uDt2p[i])
                Ycorr[i] = Ycur[i] .+ Δt/2 *(ϕ_y[i] .+ ϕ_yp[i]) .+ Δt^2 / 12 *(DvDt[i] .- DvDtp[i]) .+ Δt^3 /24 * (D2vDt2[i] .+ D2vDt2p[i])
                ϕcorr[i] = ϕcur[i] .+ Δt/2 *(DϕDt[i] .+ DϕDtp[i]) .+ Δt^2 / 12 *(D2ϕDt2[i] .- D2ϕDt2p[i]) .+ Δt^3 /24 * (D3ϕDt3[i] .+ D3ϕDt3p[i])
            end
            # Conditionally smooth the corrected surface, mirroring Dold's
            # "smooth"/"smthwd" logic: he only invokes smoothing when the
            # roughness diagnostic erm indicates the surface has picked up
            # high-wavenumber ("sawtooth") noise beyond the local error
            # tolerance - i.e. (errortol + erm)^4 > errortol - rather than
            # unconditionally every step (which is what this code did
            # previously, and which over-damps genuine short-wavelength
            # surface features). Note: this ports Dold's triggering
            # condition and his m-point binomial smoothing kernel, but not
            # his additional spatially-weighted "smthwd"/"spreadd"
            # machinery (which locally tapers the smoothing strength based
            # on nearby curvature) - that refinement is a possible future
            # addition, not implemented here.
            if smoothingval && (errortolval + erm)^4 > errortolval
                smooth!(ws.Ω_sm_temp,Xcorr,2π,N,ws.smoothCoefficients)
                smooth!(ws.Ω_sm_temp,Ycorr,0.0,N,ws.smoothCoefficients)
                smooth!(ws.Ω_sm_temp,ϕcorr,0.0,N,ws.smoothCoefficients)
            end

            # Guard against NaN/Inf (blow-up or instability) before committing
            # this step's values, mirroring Dold's per-step
            # "if (abs(erm)>5e8 .or. isnan) stop '<data is no longer
            # intelligible>'" check.
            if !(all(isfinite, Xcorr) && all(isfinite, Ycorr) && all(isfinite, ϕcorr))
                println("Data is no longer intelligible (NaN/Inf detected). Aborting simulation.")
                breaking = true
            else
                # Accept the step: always advance the true current state...
                Xcur .= Xcorr
                Ycur .= Ycorr
                ϕcur .= ϕcorr
                tcur += Δt
                prevΔt = Δt # for next step's growth-rate cap

                # ...but only append to the output series on a landing step,
                # i.e. when tcur has just reached a clean multiple of dt̃val.
                # Sub-steps taken between save points still update Xcur/
                # Ycur/ϕcur/tcur above, they just aren't recorded.
                if landingStep
                    push!(Xfull,copy(Xcur))
                    push!(Yfull,copy(Ycur))
                    push!(ϕfull,copy(ϕcur))
                    push!(t,tcur)
                    nextSaveTime += dt̃val
                end
            end
        catch e
            if e isa Union{ArgumentError,InexactError}
                println("Surface became multi-valued. Aborting simulation.")
                breaking = true
            else
                rethrow(e)
            end    
        end
    end

    # If the run ended without landing exactly on a save-time boundary (T̃val
    # need not be an exact multiple of dt̃val, and a run can also legitimately
    # reach T̃val partway through the final interval), still record the true
    # final state rather than silently dropping it from the output. Skipped
    # on an aborted (NaN/breaking) run, since Xcur/Ycur/ϕcur may not reflect
    # a valid accepted step in that case.
    if !breaking && t[end] != tcur
        push!(Xfull,copy(Xcur))
        push!(Yfull,copy(Ycur))
        push!(ϕfull,copy(ϕcur))
        push!(t,tcur)
    end

    # Reshape into matrix 
    # Note: reduce(hcat, ...) instead of hcat(Xfull...) - splatting a
    # vector with one entry per timestep as positional args is a Julia
    # anti-pattern that can dominate runtime/compile time for long runs.
    Xmatrix = permutedims(reduce(hcat, Xfull))
    Ymatrix = permutedims(reduce(hcat, Yfull))
    ϕmatrix = permutedims(reduce(hcat, ϕfull))

    # Redimensionalize variables 

    return Xmatrix*lengthScale, (Ymatrix.+MWL)*lengthScale, ϕmatrix*lengthScale*timeScale, t*timeScale
end
