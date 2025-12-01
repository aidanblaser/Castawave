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

function fixedTimeOperations(N::Int, X::AbstractVector{<:Real}, Y::AbstractVector{<:Real}, ϕ::AbstractVector{<:Real}, L::Real, h::Real)
    #=
    The fixedTimeOperations function wraps all of the necessary operations for finding the quantities needed for the next timestep. These consist of the finding the R_ξ derivative and the ϕ_ξ, ϕ_ν derivatives, from which ϕ_t is given from Bernoulli's condition. Then, the change in both R = X + iY and ϕ is known and the system can be evolved to the next timestep.
    
    It has many calls to the underlying helper functions file for clarity and compartementalization of the code.
    
    Input:
    N - number of Lagrangian particles and length of position and velocity vectors
    X - real vector of particle x-positions on the surface
    Y - real vector of particle y-positions on the surface
    ϕ - real vector of scalar velocity potential for particles on the surface

    Output:
    ϕ_x - real vector of U velocity for particles on the surface
    ϕ_y - real vector of V velocity for particles on the surface
    ϕ_D - real vector of the material derivative of ϕ from the dynamic boundary condition
    =#

    Ω = conformalMap(X .+ im*Y)
    R_ξ = DDI1(X,N,L) .+ im*DDI1(Y,N,0)
    Ω_ξ = DDI1(Ω,N,0)
    Ω_ξξ = DDI2(Ω,N,0)

    H = conformalDepth(h)

    # The matrix method described in Dold is used to find the normal derivative of the potential.
    A, B, ℵ = ABMatrices(Ω, Ω_ξ, Ω_ξξ, N, H)
    ϕ_ξ, ϕ_ν = NormalInversion(ϕ, A, ℵ, N)

    # All necessary information having been computed, we transform back to the real frame to output conveniant timestepping quantities.
    ϕ_x, ϕ_y = RealPhi(R_ξ, ϕ_ξ, ϕ_ν)
    
    # Use all this to compute up to third order time derivatives
    DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3 = TimeDerivatives(R_ξ, ϕ_x, ϕ_y, A, ℵ, N,Y)

    return ϕ_x, ϕ_y, DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3
end


function runSim(N::Int, X::AbstractVector{<:Real}, Y::AbstractVector{<:Real}, ϕ::AbstractVector{<:Real}, dt::Real, tf::Real,L::Real,h::Real; ϵ::Real = 1e-5,smoothing::Bool=true)
    #=
    The run function is the master function of the program, taking in the initial conditions for the system
    and timestepping it forward until some final time. The outputs are written into arrays
    (which can also be exported through JDL2 for larger files).
    
    Input:
    N  - number of Lagrangian particles and length of position and velocity vectors
    X  - real vector (length N) of initial particle x-positions on the surface
    Y  - real vector (length N) of initial particle y-positions on the surface
    ϕ  - real vector (length N) of initial scalar velocity potential for particles on the surface
    dt - Time step (scalar, seconds)
    tf - Final time of simulation (scalar, seconds)
    L  - Periodicity of system in physical length units (meters)
    h  - Depth of domain (scalar, meters). Defaults to 0 which means infinite depth.
    smooth - (true/false) Setting smoothing to true makes all values (X,Y,ϕ) become smoothed
              at each timestep. The smoothing is a 15 point stencil as defined in Helperfunctions.jl (function name smooth())

    Output:
    X_timeseries - real array of x-positions for particles on the surface at each time step
    X_timeseries - real array of y-positions for particles on the surface at each time step
    ϕ_timeseries - real vector of ϕ values for particles on the surface at each time step
    time         - real array of timesteps from 0 to tf with step size dt
    =#

    # Make domain 2π periodic
    L̃ = L / 2π
    hs = h / L̃


    XS = X / L̃
    YS = Y / L̃
    ϕS = ϕ / (L̃)^(3/2)

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

    while t[end] <= tf && !breaking
        try
            # smooth data if desired (and not for first timestep)
            if smoothing && length(t) > 1
                Xfull[end] = smooth(Xfull[end],N,2π)
                Yfull[end] = smooth(Yfull[end],N,0)
                ϕfull[end] = smooth(ϕfull[end],N,0)
            end

            # Determine velocities to timestep particles
            ϕ_x, ϕ_y, DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3 = fixedTimeOperations(N, Xfull[end], Yfull[end], ϕfull[end], 2π, h)




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
                Δt = (ϵ*factorial(3)/max(thirdOrderMax,fourthOrderMax))^(1/3) / 10
            else
                thirdOrderMax = max(maximum(abs.(D2uDt2)),maximum(abs.(D2vDt2)),maximum(abs.(D3ϕDt3)))
                fourthOrderMax = max(maximum(abs.(dX[:,1])),maximum(abs.(dY[:,1])),maximum(abs.(dϕ[:,1])))
                Δt = min((ϵ*factorial(3)/max(thirdOrderMax,fourthOrderMax))^(1/3),dt)
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
    return Xmatrix, Ymatrix, ϕmatrix, t
end
