#=
This file contains the main routine for initialising and running the adapted Dold flow solver. 
The program is initialized through the input parameters of the run function, which then runs the solver for 
the specified time or amount.

Package dependencies as well as the other two files needed to run the solver are initially called and included. 
The fixedTimeOperations computes all necessary quantities needed to step the model forward in time, while the 
overarching run function marches the X, Y, and ϕ Vectors in time until outputting a matrix wherein each row 
corresponds to a fixed time. 

Note that, during initialization, particle labels are taken as the indices of the X, Y, and ϕ vectors, and are 
as such evenly spaced no matter the positional values.
=#

using LinearAlgebra

include("Constants.jl")
include("Types.jl")
include("HelperFunctions.jl")




function fixedTimeOperations(X::AbstractVector{<:Real}, Y::AbstractVector{<:Real}, ϕ::AbstractVector{<:Real}, p::SimulationParameters,N::Int,H)
    #=
    The fixedTimeOperations function wraps all of the necessary operations for finding the quantities needed 
    for the next timestep. These consist of the finding the R_ξ derivative and the ϕ_ξ, ϕ_ν derivatives, 
    from which ϕ_t is given from Bernoulli's condition. Then, the change in both R = X + iY and ϕ is known and the 
    system can be evolved to the next timestep.
    
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
    (Up to 3rd order time derivatives are also computed for X, Y, ϕ)
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

    # All necessary information having been computed, we transform back to the real frame to output conveniant 
    # timestepping quantities.
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
    The fixedTimeOperations function wraps all of the necessary operations for finding the quantities needed for 
    the next timestep. These consist of the finding the R_ξ derivative and the ϕ_ξ, ϕ_ν derivatives, 
    from which ϕ_t is given from Bernoulli's condition. Then, the change in both R = X + iY and ϕ is known and 
    the system can be evolved to the next timestep.
    
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

    # All necessary information having been computed, we transform back to the real frame to output 
    # conveniant timestepping quantities.
    RealPhi!(ws, ϕ_x, ϕ_y, ϕ_ξ, ϕ_ν)

    # Compute up to third order time derivatives 
    TimeDerivatives!(ws, ℵ, ϕ_x, ϕ_y, Y, DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3)
end


function runSim(X::AbstractVector{<:Real}, Y::AbstractVector{<:Real}, ϕ::AbstractVector{<:Real}, p::SimulationParameters)
    #=
    The run function is the main function of the program, taking in the initial conditions for the system
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
             that's governed by atol (and the stability cap below); dt only
             controls how densely the returned series is sampled. Internal
             steps are taken freely between save points and simply aren't
             recorded, with the final sub-step of each interval shortened so
             the saved series lands on clean multiples of dt.
        tf - duration of simulation (seconds, will abort early if wave breaks)
        atol - absolute error tolerance (in the same nondimensional length
               units as X/Y/ϕ) on how much the predictor and corrector
               steps are allowed to disagree each internal step; this is
               what actually sets the internal timestep now (see below).
               Defaults to 1e-6. Deliberately absolute rather than
               relative: X/Y/ϕ all pass through zero somewhere on a
               periodic surface, where a relative tolerance would
               spuriously blow up right at those crossings.
        errortol - defaults to 1e-4. No longer sets the internal timestep
                   (see atol above) - only gates the conditional surface
                   smoothing below (via the roughness diagnostic erm).
        smoothing - Default to true. When enabled, gates a 15-pt smoothing
                    filter (Dold's "smthwd") that removes "sawtooth" modes
                    which tend to appear in these codes - only triggered
                    when the roughness diagnostic erm exceeds threshold,
                    and even then only applied point-by-point in proportion
                    to each point's own local roughness, so a genuinely
                    sharp but well-resolved feature isn't smoothed away
                    just because noise appeared somewhere else on the
                    surface.
        g - value of gravitational acceleration (defaults to 9.81), used
            as given (kept physical, not rescaled)

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

    # Add breaking parameter (will become true and abort code if the wave breaks)
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

    # Preallocate every scratch buffer used by the timestepper into one workspace.
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

    # --- Adaptive timestepping: predictor/corrector local error control ---
    # Δt is no longer recomputed from scratch each step via a derivative-
    # magnitude formula (Dold's sizedt) - it carries over between steps,
    # adjusted up or down by how much the predictor (Xnext/Ynext/ϕnext, a
    # 3rd-order Taylor step) and corrector (Xcorr/Ycorr/ϕcorr, a higher-
    # order trapezoidal-with-correction step) actually disagreed on the
    # last attempt. That disagreement is a direct, honest local truncation
    # error estimate - and unlike erm/dm, it's measured in the same
    # (physical-gravity, lengthScale-nondimensionalized) units everything
    # else in this solver uses, so no unit-conversion trickery is needed.
    #
    # atolVal is deliberately a pure ABSOLUTE tolerance, not relative: X,
    # Y and ϕ all legitimately pass through zero somewhere on a periodic
    # surface, and a relative tolerance would spuriously blow up right at
    # those crossings.
    atolVal = p.atol

    # Standard embedded/step-doubling step-size control (Hairer, Norsett &
    # Wanner, "Solving ODEs I" - essentially the same formula scipy's
    # solve_ivp, MATLAB's ode45, and DifferentialEquations.jl all use).
    # safety keeps the next attempt just inside tolerance rather than
    # exactly on the edge of it; facmin/facmax bound how much Δt can
    # shrink or grow in one adjustment (at most 5x either way) so it
    # doesn't oscillate step to step; pOrder is the predictor's order (3,
    # a 3rd-order Taylor step) - the exponent 1/(pOrder+1) comes from its
    # leading truncation-error term scaling like Δt^(pOrder+1).
    safety = 0.9
    facmin = 0.2
    facmax = 5.0
    pOrder = 3

    # A floor below which a step that still can't satisfy atolVal is taken
    # as the surface having genuinely broken, rather than shrinking Δt
    # forever. Deliberately chosen and adjustable - unlike the abs(erm)>
    # 5e8 threshold this replaces, which was copied from Dold's own code
    # and calibrated for his unit convention.
    minΔt = 1e-8

    # Initial guess for the very first step only - deliberately
    # conservative; the accept/reject loop below corrects it within the
    # first few steps regardless of how good or bad this guess is.
    Δt = 1e-4

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
            # Compute up to third order derivatives of X, Y, ϕ at the
            # current (already-accepted) state. Doesn't depend on Δt, so
            # it's done once per outer iteration - retries inside the
            # accuracy loop below only redo the predictor/corrector and the
            # predicted-point evaluation, not this.
            fixedTimeOperations!(ws, Xcur, Ycur, ϕcur, ϕ_x, ϕ_y,
            DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3)

            # Roughness diagnostic (Dold's "rough"/erm) - kept to gate the
            # conditional smoothing below (and, via the full per-point
            # vector ws.roughnessVec, to weight it spatially - see
            # smoothWeighted! below), not for timestep control or a hard
            # blow-up stop (both superseded by the predictor/corrector
            # error control below, which measures accuracy directly
            # instead of inferring it from a raw derivative magnitude).
            erm = roughnessVector!(ws.roughnessVec, D2uDt2, D2vDt2, D3ϕDt3, N, ws.roughnessCoefficients)

            # Numerical-stability threshold on the timestep, ported from Dold's
            # "timstp"/strong-instability check (Dold 1992, J. Comp. Phys. 103,
            # 90-115). Explicit surface-tracking schemes like this one become
            # strongly unstable above a step size set by the local balance of
            # surface curvature/acceleration against gravity - independent of
            # the accuracy-based error control below, which can look fine
            # while a run is quietly heading into instability (especially as
            # a wave steepens toward breaking). stabilityFactorVal defaults
            # to 0.6, matching Dold's own recommended (most conservative)
            # setting. Applied to Δt (the carried-forward, "natural" step
            # size) once per outer iteration, before any retries.
            stabilityMax = 0.0
            @inbounds for i ∈ 1:N
                tp = abs(ws.X_ξ[i]*(DvDt[i]+gravity) - ws.Y_ξ[i]*DuDt[i]) / (ws.X_ξ[i]^2 + ws.Y_ξ[i]^2)
                stabilityMax = max(stabilityMax, tp)
            end
            Δt = min(Δt, stabilityFactorVal / sqrt(stabilityMax))

            # --- Accuracy-controlled predictor/corrector retry loop ---
            # Δt is the carried-forward "natural" step size, evolved only
            # by the accuracy-based growth/shrink factor below - it's never
            # overwritten by the landing-clip logic. stepΔt is the (possibly
            # landing-clipped) value actually used for this attempt; it's
            # recomputed from the current Δt every retry, so a step that
            # gets rejected and shrunk also gets its landing-clip check
            # redone against the smaller Δt.
            accepted = false
            landingStep = false
            stepΔt = Δt
            while !accepted
                if Δt < minΔt
                    println("Data is no longer intelligible (Δt fell below $(minΔt) while still failing the position/potential error tolerance atol=$(atolVal)). Aborting simulation.")
                    breaking = true
                    break
                end

                # Clip the trial step so it approaches the next save-time
                # (dt̃val-spaced) boundary, exactly as before: Dold's own
                # "pts" printout interval - a fixed interval, no breaking-
                # triggered densification - including his nicedt "homing in
                # smoothly from as many as almost 4 steps distant" (the
                # step size eases down over the last few steps rather than
                # clipping one abrupt final step). This is purely local to
                # this attempt - it never feeds back into Δt itself.
                stepΔt = Δt
                landingStep = false
                if nextSaveTime <= T̃val
                    if tcur + 1.1*stepΔt >= nextSaveTime
                        stepΔt = nextSaveTime - tcur
                        landingStep = true
                    elseif tcur + 1.9*stepΔt >= nextSaveTime
                        stepΔt = 0.5*(nextSaveTime - tcur)
                    elseif tcur + 2.9*stepΔt >= nextSaveTime
                        stepΔt = (1.0/3.0)*(nextSaveTime - tcur)
                    elseif tcur + 3.9*stepΔt >= nextSaveTime
                        stepΔt = 0.25*(nextSaveTime - tcur)
                    end
                end

                # Predictor: 3rd-order Taylor step.
                @inbounds for i ∈ 1:N
                    Xnext[i] = Xcur[i] + ϕ_x[i]*stepΔt + stepΔt^2/2*DuDt[i] + stepΔt^3/6*D2uDt2[i]
                    Ynext[i] = Ycur[i] + stepΔt*ϕ_y[i] + stepΔt^2/2*DvDt[i] + stepΔt^3/6*D2vDt2[i]
                    ϕnext[i] = ϕcur[i] + stepΔt*DϕDt[i] + stepΔt^2/2*D2ϕDt2[i] + stepΔt^3/6*D3ϕDt3[i]
                end

                # Estimate derivatives at the predicted surface
                fixedTimeOperations!(ws, Xnext, Ynext, ϕnext, ϕ_xp, ϕ_yp,
                DϕDtp, DuDtp, DvDtp, D2ϕDt2p, D2uDt2p, D2vDt2p, D3ϕDt3p)

                # Corrector: trapezoidal average of derivatives at the
                # predicted and current surfaces, with Taylor correction
                # terms - higher order than the predictor, so the two
                # together form an embedded pair.
                @inbounds for i ∈ 1:N
                    Xcorr[i] = Xcur[i] + stepΔt/2*(ϕ_x[i]+ϕ_xp[i]) + stepΔt^2/12*(DuDt[i]-DuDtp[i]) + stepΔt^3/24*(D2uDt2[i]+D2uDt2p[i])
                    Ycorr[i] = Ycur[i] + stepΔt/2*(ϕ_y[i]+ϕ_yp[i]) + stepΔt^2/12*(DvDt[i]-DvDtp[i]) + stepΔt^3/24*(D2vDt2[i]+D2vDt2p[i])
                    ϕcorr[i] = ϕcur[i] + stepΔt/2*(DϕDt[i]+DϕDtp[i]) + stepΔt^2/12*(D2ϕDt2[i]-D2ϕDt2p[i]) + stepΔt^3/24*(D3ϕDt3[i]+D3ϕDt3p[i])
                end

                # Local error estimate: how much the predictor and
                # corrector disagree, in absolute (X,Y,ϕ) units, relative
                # to atolVal - computed BEFORE any smoothing is applied
                # below, so smoothing (which deliberately damps
                # high-wavenumber noise) can't mask a genuinely inaccurate
                # step.
                errX = 0.0; errY = 0.0; errϕ = 0.0
                @inbounds for i ∈ 1:N
                    errX = max(errX, abs(Xcorr[i]-Xnext[i]))
                    errY = max(errY, abs(Ycorr[i]-Ynext[i]))
                    errϕ = max(errϕ, abs(ϕcorr[i]-ϕnext[i]))
                end
                errRatio = max(errX, errY, errϕ) / atolVal

                # Standard embedded/step-doubling step-size update (Hairer,
                # Norsett & Wanner) - used both to shrink Δt on rejection
                # and to grow/shrink it modestly on acceptance.
                factor = clamp(safety*errRatio^(-1/(pOrder+1)), facmin, facmax)
                Δt = Δt * factor

                if errRatio <= 1.0
                    accepted = true
                end
                # else: rejected - loop back and retry with the smaller Δt
                # (stepΔt, and its landing-clip check, are both recomputed
                # from it at the top of the loop).
            end

            if breaking
                break
            end

            # Conditionally smooth the corrected surface, mirroring Dold's
            # "smthwd" logic: he only invokes smoothing when the roughness
            # diagnostic erm indicates the surface has picked up
            # high-wavenumber ("sawtooth") noise beyond the local error
            # tolerance - i.e. (errortol + erm)^4 > errortol - rather than
            # unconditionally every step. And, per Dold's own documentation
            # ("smoothing is applied only if the effect of roughness
            # exceeds |precision|^(1/4), and then only where the roughness
            # is most significant"), the correction is also weighted
            # per-point by roughnessWeight! before being applied, so a
            # genuinely sharp but well-resolved feature (e.g. a plunging
            # jet tip) isn't flattened just because smoothing fired
            # elsewhere on the surface - only points whose own local
            # roughness signal is high get corrected at full strength.
            if smoothingval && (errortolval + erm)^4 > errortolval
                roughnessWeight!(ws.roughnessWeight, ws.roughnessWeightTemp, ws.roughnessVec, erm, N, ws.spreadCoefficients)
                smoothWeighted!(ws.Ω_sm_temp,Xcorr,2π,N,ws.smoothCoefficients,ws.roughnessWeight)
                smoothWeighted!(ws.Ω_sm_temp,Ycorr,0.0,N,ws.smoothCoefficients,ws.roughnessWeight)
                smoothWeighted!(ws.Ω_sm_temp,ϕcorr,0.0,N,ws.smoothCoefficients,ws.roughnessWeight)
            end

            # Guard against NaN/Inf (blow-up or instability) before committing
            # this step's values. Mostly a redundant final safety net now
            # (a NaN in Xcorr/Xnext would typically make errRatio itself
            # NaN, and NaN <= 1.0 is false in Julia, so the retry loop above
            # would already reject and shrink toward minΔt on its own) but
            # kept as a direct, unambiguous check rather than relying on
            # that indirectly.
            if !(all(isfinite, Xcorr) && all(isfinite, Ycorr) && all(isfinite, ϕcorr))
                println("Data is no longer intelligible (NaN/Inf detected). Aborting simulation.")
                breaking = true
            else
                # Accept the step: always advance the true current state...
                Xcur .= Xcorr
                Ycur .= Ycorr
                ϕcur .= ϕcorr
                tcur += stepΔt

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

    println(tcur)

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
    Xmatrix = permutedims(reduce(hcat, Xfull))
    Ymatrix = permutedims(reduce(hcat, Yfull))
    ϕmatrix = permutedims(reduce(hcat, ϕfull))

    # Redimensionalize variables 

    return Xmatrix*lengthScale, (Ymatrix.+MWL)*lengthScale, ϕmatrix*lengthScale*timeScale, t*timeScale
end
