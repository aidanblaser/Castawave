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
        return new(L,h,dt,T,lengthScale,timeScale,h̃,dt̃,T̃,errortol,smoothing,g)
    end 
end


function fixedTimeOperations(X::AbstractVector{<:Real}, Y::AbstractVector{<:Real}, ϕ::AbstractVector{<:Real}, p::SimulationParameters)
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
    N = length(X)
    h̃ = p.h̃

    Ω = conformalMap(X .+ im*Y)
    # Requires to be non-dimensionsionalized
    R_ξ = DDI1(X,2π) .+ im*DDI1(Y,0)
    Ω_ξ = DDI1(Ω,0)
    Ω_ξξ = DDI2(Ω,0)

    H = conformalDepth(h)

    # The matrix method described in Dold is used to find the normal derivative of the potential.
    A, B, ℵ = ABMatrices(Ω, Ω_ξ, Ω_ξξ, H)
    ϕ_ξ, ϕ_ν = NormalInversion(ϕ, A, ℵ)

    # All necessary information having been computed, we transform back to the real frame to output conveniant timestepping quantities.
    ϕ_x, ϕ_y = RealPhi(R_ξ, ϕ_ξ, ϕ_ν)
    
    # Use all this to compute up to third order time derivatives
    DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3 = TimeDerivatives(R_ξ, ϕ_x, ϕ_y, A, ℵ,Y,p)

    return ϕ_x, ϕ_y, DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3
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
    hS = h / lengthScale

    # Add breaking parameter 
    breaking = false

    #MWL
    #MWL =  sum(DDI1(X,N,L,1).*Y)/L

    # Initialize time vector
    t = [0.0]

    # Create and initialize the timeseries fields. 
    Xfull = Vector{Vector{Float64}}()
    Yfull = Vector{Vector{Float64}}()
    ϕfull = Vector{Vector{Float64}}()
    push!(Xfull,XS)
    push!(Yfull,YS)
    push!(ϕfull,ϕS)

    # Preallocate derivative matrices
    dX = Array{Float64}(undef,N,5)
    dY = Array{Float64}(undef,N,5)
    dϕ = Array{Float64}(undef,N,5)

    while t[end] <= p.T̃ && !breaking
        try
            # smooth data if desired (and not for first timestep)
            if p.smoothing && length(t) > 1
                Xfull[end] = smooth(Xfull[end],2π)
                Yfull[end] = smooth(Yfull[end],0)
                ϕfull[end] = smooth(ϕfull[end],0)
            end

            # Determine velocities to timestep particles
            ϕ_x, ϕ_y, DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3 = fixedTimeOperations(Xfull[end], Yfull[end], ϕfull[end],p)

            # For each point
            Xnext = similar(XS)
            Ynext = similar(YS)
            ϕnext = similar(ϕS)
            for i ∈ 1:N
                # From third order derivatives, extrapolate higher order using Lagrange polynomials
                # Handle first few timesteps separately
                if length(t) < 5
                    dX[i,:] = LagrangeInterpolant([Xfull[k][i] for k in 1:length(t)],t,length(t)-1)
                    dY[i,:] = LagrangeInterpolant([Yfull[k][i] for k in 1:length(t)],t,length(t)-1)
                    dϕ[i,:] = LagrangeInterpolant([ϕfull[k][i] for k in 1:length(t)],t,length(t)-1)
                else
                    dX[i,:] = LagrangeInterpolant([Xfull[k][i] for k in length(t)-4:length(t)],t[end-4:end],5)
                    dY[i,:] = LagrangeInterpolant([Yfull[k][i] for k in length(t)-4:length(t)],t[end-4:end],5)
                    dϕ[i,:] = LagrangeInterpolant([ϕfull[k][i] for k in length(t)-4:length(t)],t[end-4:end],5)
                end
            end

            # Determing Timestep
            if length(t) < 5
                thirdOrderMax = max(maximum(abs.(D2uDt2)),maximum(abs.(D2vDt2)),maximum(abs.(D3ϕDt3)))
                fourthOrderMax = max(maximum(abs.(dX[:,1])),maximum(abs.(dY[:,1])),maximum(abs.(dϕ[:,1])))
                # Take smaller timesteps initially
                Δt = (p.errortol*factorial(3)/max(thirdOrderMax,fourthOrderMax))^(1/3) / 10
            else
                thirdOrderMax = max(maximum(abs.(D2uDt2)),maximum(abs.(D2vDt2)),maximum(abs.(D3ϕDt3)))
                fourthOrderMax = max(maximum(abs.(dX[:,1])),maximum(abs.(dY[:,1])),maximum(abs.(dϕ[:,1])))
                Δt = min((p.errortol*factorial(3)/max(thirdOrderMax,fourthOrderMax))^(1/3),p.dt̃)
                # Minimum timestep 1e-5
                Δt = max(Δt,1e-5)
            end

            for i ∈ 1:N

                # Use all these derivatives to get next timestep
                Xnext[i] = Xfull[end][i] +Δt * (ϕ_x[i]) +Δt^2/factorial(2)*DuDt[i] + Δt^3/factorial(3)*D2uDt2[i] +
                sum([Δt^n/factorial(n)*dX[n-3] for n=4:8])
                Ynext[i] = Yfull[end][i] +Δt * (ϕ_y[i]) + Δt^2/factorial(2)*DvDt[i] +Δt^3/factorial(3)*D2vDt2[i] +
                sum([Δt^n/factorial(n)*dY[n-3] for n=4:8])
                ϕnext[i] = ϕfull[end][i] +Δt * (DϕDt[i]) +Δt^2/factorial(2)*D2ϕDt2[i] +Δt^3/factorial(3)*D3ϕDt3[i]+
                sum([Δt^n/factorial(n)*dϕ[n-3] for n=4:8])
            end
            # Append these values to the result
            push!(Xfull,Xnext)
            push!(Yfull,Ynext)
            push!(ϕfull,ϕnext)
            push!(t,t[end] + Δt)
        catch e
            if e isa ArgumentError 
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

    return Xmatrix*lengthScale, Ymatrix*lengthScale, ϕmatrix*lengthScale*timeScale, t*timeScale
end
