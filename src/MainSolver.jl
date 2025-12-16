#=
This file contains the main routine for initialising and running the adapted Dold flow solver. The program is initialized through the input parameters of the run function, which then runs the solver for the specified time or amount.

Package dependencies as well as the other two files needed to run the solver are initially called and included. The fixedTimeOperations computes all necessary quantities needed to step the model forward in time, while the overarching run function marches the X, Y, and ϕ Vectors in time until outputting a matrix wherein each row corresponds to a fixed time. 

Note that, during initialization, particle labels are taken as the indices of the X, Y, and ϕ vectors, and are as such evenly spaced no matter the positional values.
=#

using DrWatson
@quickactivate "Castawave"
using LinearAlgebra

include(projectdir()*"/src/Constants.jl")
include(projectdir()*"/src/HelperFunctions.jl")

#= TODO 
- Better implement time stepping so it lands on intervals of dt (even if it takes steps between)
- Better implement code aborting condition (maybe based on min timestep?)
=#

# Create a structure for the parameters so it knows it's invariant
struct SimulationParameters

    # Physical parameters
    L::Float64 # Length of domain in meters (spatial periodicity)
    h::Float64 # Depth of fluid at rest in meters
    dt::Float64 # Maximum timestep in seconds 
    T::Float64 # Desired duration of simulation in seconds (may abort early if wave breaks)

    # Parameters re-scaled to have L=2π
    lengthScale::Float64
    timeScale::Float64
    h̃::Float64 
    dt̃::Float64 
    T̃::Float64

    # Parameters with default values 
    errortol::Float64
    smoothing::Bool 
    g::Float64
end

function SimulationParameters(L,h,dt,T;errortol=1e-5,smoothing=true,g=9.81)
        # Convert variables to explicit types 
        L = Float64(L)
        h = Float64(h)
        dt = Float64(dt)
        T = Float64(T)
        errortol = Float64(errortol)
        smoothing = Bool(smoothing)
        g = Float64(g)

        # Parameters scaled to have 2π periodicity
        lengthScale = Float64(L/2π)
        timeScale = sqrt(lengthScale)
        h̃ = h/lengthScale
        dt̃ = dt/timeScale
        T̃ = T/timeScale
        return SimulationParameters(L, h, dt, T,
                               lengthScale, timeScale, h̃, dt̃, T̃,
                               errortol, smoothing, g)
end 



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

function fixedTimeOperations!(Ω,X,Y,X_ξ,N,c1,c2,Y_ξ,Ω_ξ,Ω_ξξ,inverseds2,A,B,C,Cinter,ΔΩ,
    H,ϕ_ξ,ϕ_ξξ,ϕ_ν,b,ϕ,ϕ_x,ϕ_y,
    DϕDt,DuDt,DvDt,D2ϕDt2,D2uDt2,D2vDt2,D3ϕDt3,
    ϕ_t,ϕ_tξ,ϕ_tξξ,ϕ_tν,u_t,u_tξ,v_t,v_tξ,u_ξ,v_ξ,u_x,v_x,u_tx,v_tx,u_xξ,v_xξ,u_xx,v_xx,
    ϕ_tt,ϕ_ttξ,ϕ_ttξξ,ϕ_ttν,u_tt,v_tt,g)
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
    DDI1!(X_ξ,X,2π,N,c1)
    DDI1!(Y_ξ,Y,0.0,N,c1)
    DDI1!(Ω_ξ,Ω,zero(eltype(Ω)),N,c1)
    DDI2!(Ω_ξξ,Ω,zero(eltype(Ω)),N,c2)

    # Computing inverse spacing (used a lot)
    @inbounds for i ∈ 1:N 
        inverseds2[i] = (X_ξ[i]^2 + Y_ξ[i]^2)^(-1)
    end 

    # The matrix method described in Dold is used to find the normal derivative of the potential.
    ℵ = ABMatrices!(A,B,C,Cinter,ΔΩ, Ω, Ω_ξ, Ω_ξξ, H,N)
    NormalInversion!(ϕ_ξ,ϕ_ξξ,ϕ_ν,b,ϕ, A, ℵ,N,c1,c2)

    # All necessary information having been computed, we transform back to the real frame to output conveniant timestepping quantities.
    RealPhi!(ϕ_x,ϕ_y,X_ξ,Y_ξ,inverseds2,ϕ_ξ, ϕ_ν,N)

    # Checked up to here 

    # Compute up to third order time derivatives 
    TimeDerivatives!(DϕDt,DuDt,DvDt,D2ϕDt2,D2uDt2,D2vDt2,D3ϕDt3,
    X_ξ,Y_ξ,ϕ_x,ϕ_y,A,b,ℵ,Y,N,c1,c2,
    ϕ_t,ϕ_tξ,ϕ_tξξ,ϕ_tν,u_t,u_tξ,v_t,v_tξ,u_ξ,v_ξ,u_x,v_x,u_tx,v_tx,u_xξ,v_xξ,u_xx,v_xx,
    ϕ_tt,ϕ_ttξ,ϕ_ttξξ,ϕ_ttν,u_tt,v_tt,inverseds2,g)
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
        dt - desired timestep (seconds)
        tf - duration of simulation (seconds, will abort early if wave breaks)
        errortol - error tolerance, defaults to 1e-5. Lower values mean greater accuracy. Reduces timestep based on nonlinearity
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
    
    # Create and initialize the timeseries fields. 
    Xfull = Vector{Vector{Float64}}()
    Yfull = Vector{Vector{Float64}}()
    ϕfull = Vector{Vector{Float64}}()
    push!(Xfull,XS)
    push!(Yfull,YS .- MWL)
    push!(ϕfull,ϕS)

    c1 =[2100.0, -600.0, 150.0, -25.0, 2.0] ./ 2520.0;


    c2 = [-9220.75, 5250.0, -750.0, 125.0, -15.625, 1.0]./3150;

    Smoothcoefficients =  [3432.0, -3003.0, 2002.0, -1001.0, 364.0, -91.0, 14.0, -1.0]./(2^14);

    # Preallocate all intermediate variables 
    Ω_sm_temp = similar(X)
    X_ξ = similar(X)
    Y_ξ = similar(X)
    Xnext = similar(X)
    Ynext = similar(X)
    ϕnext = similar(X)
    Xcorr = similar(X)
    Ycorr = similar(X)
    ϕcorr = similar(X)
    ϕ_x = similar(X)
    ϕ_y = similar(X)
    DϕDt = similar(X)
    DuDt= similar(X)
    DvDt= similar(X)
    D2ϕDt2 = similar(X)
    D2uDt2 = similar(X)
    D2vDt2 = similar(X)
    D3ϕDt3 = similar(X)
    ϕ_xp = similar(X)
    ϕ_yp = similar(X)
    DϕDtp = similar(X)
    DuDtp = similar(X)
    DvDtp = similar(X)
    D2ϕDt2p = similar(X)
    D2uDt2p = similar(X)
    D2vDt2p = similar(X)
    D3ϕDt3p = similar(X)
    ϕ_tξξ = similar(X)
    ϕ_ttξξ = similar(X)
    # Temporary variables 
    ϕ_ξ = similar(X)
    ϕ_ξξ = similar(X)
    ϕ_ν = similar(X)
    b = similar(X)
    inverseds2 = similar(X)
    u_t = similar(X)
    u_tξ = similar(X)
    v_t = similar(X)
    v_tξ = similar(X)
    ϕ_t = similar(X)
    ϕ_tξ = similar(X)
    ϕ_tν = similar(X)
    u_ξ = similar(X)
    v_ξ = similar(X)
    u_x = similar(X)
    v_x = similar(X)
    u_tx = similar(X)
    v_tx = similar(X)
    u_xξ = similar(X)
    v_xξ = similar(X)
    u_xx = similar(X)
    v_xx = similar(X)
    ϕ_tt  = similar(X)
    ϕ_ttξ  = similar(X)
    ϕ_ttν  = similar(X)
    u_tt  = similar(X)
    v_tt  = similar(X)
    # Preallocate all complex variables 
    Ω = Vector{ComplexF64}(undef,N)
    Ω_ξ = similar(Ω)
    Ω_ξξ = similar(Ω)
    R_ξ = similar(Ω)
    A = Matrix{Float64}(undef,N,N)
    B = similar(A)
    Cinter = similar(A)
    C = Matrix{ComplexF64}(undef,N,N)
    ΔΩ = similar(C)

    T̃val = p.T̃
    smoothingval = p.smoothing 
    errortolval = p.errortol
    dt̃val = p.dt̃


    while t[end] <= T̃val && !breaking
        try
            # smooth data if desired (and not for first timestep)
            if smoothingval && length(t) > 1
                smooth!(Ω_sm_temp,Xfull[end],2π,N,Smoothcoefficients)
                smooth!(Ω_sm_temp,Yfull[end],0.0,N,Smoothcoefficients)
                smooth!(Ω_sm_temp,ϕfull[end],0.0,N,Smoothcoefficients)
            end

            # Compute up to third order derivatives of X, Y, ϕ 
            fixedTimeOperations!(Ω,Xfull[end],Yfull[end],X_ξ,N,c1,c2,Y_ξ,Ω_ξ,Ω_ξξ,inverseds2,A,B,C,Cinter,ΔΩ,
            H,ϕ_ξ,ϕ_ξξ,ϕ_ν,b,ϕfull[end],ϕ_x,ϕ_y,
            DϕDt,DuDt,DvDt,D2ϕDt2,D2uDt2,D2vDt2,D3ϕDt3,
            ϕ_t,ϕ_tξ,ϕ_tξξ,ϕ_tν,u_t,u_tξ,v_t,v_tξ,u_ξ,v_ξ,u_x,v_x,u_tx,v_tx,u_xξ,v_xξ,u_xx,v_xx,
            ϕ_tt,ϕ_ttξ,ϕ_ttξξ,ϕ_ttν,u_tt,v_tt,gravity)


            # Determing Adaptive Timestep
            thirdOrderMax = max(maximum(abs.(D2uDt2)),maximum(abs.(D2vDt2)),maximum(abs.(D3ϕDt3)))
            Δt = min((errortolval*factorial(3)/thirdOrderMax)^(1/3),dt̃val)
            # Minimum timestep 1e-4
            Δt = max(Δt,1e-4)

            
            # Use derivatives up to third order to make a predictor step 
            @inbounds for i ∈ 1:N
                Xnext[i] = Xfull[end][i] .+ ϕ_x[i] * Δt .+ Δt^2/2*DuDt[i] .+ Δt^3/6*D2uDt2[i] 
                Ynext[i] = Yfull[end][i] .+ Δt * (ϕ_y[i]) .+ Δt^2/2*DvDt[i] .+Δt^3/6*D2vDt2[i]
                ϕnext[i] = ϕfull[end][i] .+ Δt * (DϕDt[i]) .+ Δt^2/2*D2ϕDt2[i] .+Δt^3/6*D3ϕDt3[i]
            end

            # Estimate derivatives at the predicted surface
            fixedTimeOperations!(Ω,Xnext,Ynext,X_ξ,N,c1,c2,Y_ξ,Ω_ξ,Ω_ξξ,inverseds2,A,B,C,Cinter,ΔΩ,
            H,ϕ_ξ,ϕ_ξξ,ϕ_ν,b,ϕnext,ϕ_xp,ϕ_yp,
            DϕDtp,DuDtp,DvDtp,D2ϕDt2p,D2uDt2p,D2vDt2p,D3ϕDt3p,
            ϕ_t,ϕ_tξ,ϕ_tξξ,ϕ_tν,u_t,u_tξ,v_t,v_tξ,u_ξ,v_ξ,u_x,v_x,u_tx,v_tx,u_xξ,v_xξ,u_xx,v_xx,
            ϕ_tt,ϕ_ttξ,ϕ_ttξξ,ϕ_ttν,u_tt,v_tt,gravity)

            # Use predictor-corrector to average derivatives at predicted surface and current surface (like trapezoidal rule)
            @inbounds for i ∈ 1:N
                Xcorr[i] = Xfull[end][i] .+ Δt/2 *(ϕ_x[i] .+ ϕ_xp[i]) .+ Δt^2 / 12 *(DuDt[i] .- DuDtp[i]) .+ Δt^3 /24 * (D2uDt2[i] .+ D2uDt2p[i])
                Ycorr[i] = Yfull[end][i] .+ Δt/2 *(ϕ_y[i] .+ ϕ_yp[i]) .+ Δt^2 / 12 *(DvDt[i] .- DvDtp[i]) .+ Δt^3 /24 * (D2vDt2[i] .+ D2vDt2p[i])
                ϕcorr[i] = ϕfull[end][i] .+ Δt/2 *(DϕDt[i] .+ DϕDtp[i]) .+ Δt^2 / 12 *(D2ϕDt2[i] .- D2ϕDt2p[i]) .+ Δt^3 /24 * (D3ϕDt3[i] .+ D3ϕDt3p[i])
            end
            # Append these values to the result
            push!(Xfull,copy(Xcorr))
            push!(Yfull,copy(Ycorr))
            push!(ϕfull,copy(ϕcorr))
            push!(t,t[end] + Δt)
        catch e
            if e isa Union{ArgumentError,InexactError}
                println("Surface became multi-valued. Aborting simulation.")
                breaking = true
            else
                rethrow(e)
            end    
        end
    end
    # Reshape into matrix 
    Xmatrix = permutedims(hcat(Xfull...))
    Ymatrix = permutedims(hcat(Yfull...))
    ϕmatrix = permutedims(hcat(ϕfull...))

    # Redimensionalize variables 

    return Xmatrix*lengthScale, (Ymatrix.+MWL)*lengthScale, ϕmatrix*lengthScale*timeScale, t*timeScale
end
