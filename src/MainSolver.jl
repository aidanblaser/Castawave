#=
This file contains the main routine for initialising and running the adapted Dold flow solver. The program is initialized through the input parameters of the run function, which then runs the solver for the specified time or amount.

Package dependencies as well as the other two files needed to run the solver are initially called and included. The fixedTimeOperations computes all necessary quantities needed to step the model forward in time, while the overarching run function marches the X, Y, and ϕ Vectors in time until outputting a matrix wherein each row corresponds to a fixed time. 

Note that, during initialization, particle labels are taken as the indices of the X, Y, and ϕ vectors, and are as such evenly spaced no matter the positional values.
=#

using DrWatson
@quickactivate "Castawave"
using LinearAlgebra 
using JLD2

include(projectdir()*"/src/Constants.jl")
include(projectdir()*"/src/HelperFunctions.jl")

function fixedTimeOperations(N::Int, X::Vector, Y::Vector, ϕ::Vector, L, h::Float64, smooth=false)
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

    R = [X[i] + im * Y[i] for i in 1:N]

    Ω, r, θ = conformalMap(R)
    

    H = conformalDepth(h)

    # Ω = smooth(N, Ω)
    # ϕ = smooth(N, ϕ, 1)

    # The necessary derivatives are calculated from the 11-point interpolation in the conformal frame, making them implicitly periodic and removing the need to offset the x-domain length.
    Ω_ξ = DDI1(Ω, N)
    Ω_ξξ = DDI2(Ω, N)
    R_ξ = (im ./ Ω) .* Ω_ξ

    # The matrix method described in Dold is used to find the normal derivative of the potential.
    A, B, ℵ = ABMatrices(Ω, Ω_ξ, Ω_ξξ, N, H)
    ϕ_ξ, ϕ_ν = NormalInversion(ϕ, A, ℵ, N)

    # All necessary information having been computed, we transform back to the real frame to output conveniant timestepping quantities.
    ϕ_D, ϕ_t = PhiTimeDer(R_ξ, ϕ_ξ, ϕ_ν, Y, L)
    ϕ_x, ϕ_y = RealPhi(R_ξ, ϕ_ξ, ϕ_ν)

    # The quantities for second-order Taylor series timestepping are the material derivatives. Not needed for the RK4 scheme.
    # U_D, V_D, ϕ_DD = PhiSecondTimeDer(R_ξ, ϕ_x, ϕ_y, ϕ_t, A, ℵ, N, ϕ_ξ, ϕ_ν)

    return ϕ_x, ϕ_y, ϕ_D
end


function run(N::Int, X::Vector, Y::Vector, ϕ::Vector, dt::Float64, tf::Float64,L = 2π,h=0.0,smooth = false)
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

    

    t = 0.0 #initial time
    i = 1   # timeseries index counter.
    l = ceil(Int, tf/dt) + 2   # Number of timesteps

    # Create and initialize the timeseries fields. 
    X_timeseries = zeros((l, N))
    Y_timeseries = zeros((l, N))
    ϕ_timeseries = zeros((l, N))
    time = collect(0:dt:tf)

    X_timeseries[1,:] = XS
    Y_timeseries[1,:] = YS
    ϕ_timeseries[1,:] = ϕS

    # Running the system until final time is reached. Modify function to take parameters for whether or not mean water level should be computed, and which time-scheme to use.
    while t < tf
        # for j in 1:N
        #     if X_timeseries[i, j] > 2*pi - 0.1
        #         X_timeseries[i, j] -= 2*pi
        #     end
        # end

        if mod(i, 5) == 4
            ϕ_x, ϕ_y, ϕ_D, wl[i], ta[i] = fixedTimeOperations(N, X_timeseries[i,:], Y_timeseries[i,:], ϕ_timeseries[i,:], L, hs)
        else
            ϕ_x, ϕ_y, ϕ_D, wl[i], ta[i] = fixedTimeOperations(N, X_timeseries[i,:], Y_timeseries[i,:], ϕ_timeseries[i,:], L, hs)
        end

        X_timeseries[i+1,:], Y_timeseries[i+1,:], ϕ_timeseries[i+1,:] = RK4i(dt/sqrt(L̃), fixedTimeOperations, N, X_timeseries[i,:], Y_timeseries[i,:], ϕ_timeseries[i,:], L, hs)

        i += 1
        t += dt
    end

    return X_timeseries.* L̃, Y_timeseries .*L̃, ϕ_timeseries .* L̃^(3/2),time
end
