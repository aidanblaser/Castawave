#=
This file contains the main routine for initialising and running the adapted Dold flow solver. The program is initialized through the input parameters of the run function, which then runs the solver for the specified time or amount.

Package dependencies as well as the other two files needed to run the solver are initially called and included. The fixedTimeOperations computes all necessary quantities needed to step the model forward in time, while the overarching run function marches the X, Y, and ϕ Vectors in time until outputting a matrix wherein each row corresponds to a fixed time. 

Note that, during initialization, particle labels are taken as the indices of the X, Y, and ϕ vectors, and are as such evenly spaced no matter the positional values.
=#

using DrWatson
@quickactivate "Castawave"
using LinearAlgebra 
using JLD2
using DifferentialEquations

include(projectdir()*"/src/Constants.jl")
include(projectdir()*"/src/HelperFunctions.jl")

function fixedTimeOperations(N::Int, X::Vector, Y::Vector, ϕ::Vector, L::Float64, h::Float64,smoothing=false)
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
    R_ξ = DDI1(X,N,L,1) .+ im*DDI1(Y,N,0,1)
    R_ξξ = DDI2(X,N,L,1) .+ im*DDI2(Y,N,0,1)
    Ω_ξ = -im*R_ξ.*Ω
    Ω_ξξ = -Ω.*R_ξ.^2 .- im*Ω.*R_ξξ

    H = conformalDepth(h)
    

    # The necessary derivatives are calculated from the 11-point interpolation in the conformal frame, making them implicitly periodic and removing the need to offset the x-domain length.
    # Ω_ξ = DDI1(Ω, N)
    # Ω_ξξ = DDI2(Ω, N)
    # R_ξ = (im .* Ω_ξ) ./ Ω

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

function TimeStep(du,u,p,t)
    #=
    Timestepping ODE to be used for DifferentialEquations julia package

    Input:
    du = (dX,dY,dϕ) - 3N length vector showing the evolution of each point's X,Y, and ϕ values
    u = (X,Y,ϕ) - 3N length vector showing values of coordinate
    p = (N, L, h, smoothing) - vector of parameters
    t = tf - final time to simulate
    =#
    
    N = Int(p[1]) # number of points
    u1 = u[1:N];
    u2 = u[N+1:2*N];
    u3 = u[2*N+1:3*N];
    L = p[2]
    h = p[3]
    smoothing = Bool(p[4])
    diffu1, diffu2, diffu3 = fixedTimeOperations(N,u1,u2,u3,L,h,smoothing)
    du[1:N] = diffu1
    du[N+1:2*N] = diffu2
    du[2*N+1:3*N] = diffu3
end


function runSim(N::Int, X::Vector, Y::Vector, ϕ::Vector, dt::Float64, tf::Float64,L = 2π,h=0.0;smoothing=false,MWL_reset=false,alg = Vern9(),reltol=1e-6)
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

    
    #MWL
    MWL =  sum(DDI1(X,N,L,1).*Y)/L

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

    initial = vcat(XS,YS,ϕS)
    params = [N,L,h,smoothing,MWL]

    # Callback 1: smoothing 
    function affect!(integrator)
        N = Int(integrator.p[1])
        L = integrator.p[2]
        integrator.u[1:N] = smooth(N,integrator.u[1:N],L,1)
        integrator.u[N+1:2*N] = smooth(N,integrator.u[N+1:2*N],0,1)
        integrator.u[2*N+1:3*N] = smooth(N,integrator.u[2*N+1:3*N],0,1)
    end
    condition(u,t,integrator) = smoothing
    cb1 = DiscreteCallback(condition, affect!)
    # Callback 2: Resetting MWL 
    function affectMWL!(integrator)
        N = Int(integrator.p[1])
        L = integrator.p[2]
        MWL = integrator.p[5]
        Xξ = DDI1(integrator.u[1:N],N,L,1)
        # Reset MWL 
        integrator.u[N+1:2*N] .-= (sum(Xξ.*integrator[N+1:2*N])/(L) - MWL)
    end
    condition_MWL(u,t,integrator) = MWL_reset
    cb2 = DiscreteCallback(condition_MWL,affectMWL!)

    # Define array of callbacks
    cb = CallbackSet(cb1,cb2)

    prob = ODEProblem(TimeStep,initial,tf/sqrt(L̃),params)
    #sol = solve(prob,alg,reltol=1e-8,abstol=1e-8,dtmin=1e-5,callback=cb)
    sol = solve(prob,alg,reltol=reltol,abstol=1e-8,dt=dt,callback=cb)

    # Running the system until final time is reached. Modify function to take parameters for whether or not mean water level should be computed, and which time-scheme to use.
    # while t < tf
    #    prob = ODEProblem(fixedTimeOperations,)

    #     X_timeseries[i+1,:], Y_timeseries[i+1,:], ϕ_timeseries[i+1,:] = RK4i(dt/sqrt(L̃), fixedTimeOperations, N, X_timeseries[i,:], Y_timeseries[i,:], ϕ_timeseries[i,:], L, hs,false)

    #     i += 1
    #     t += dt
    # end

    return sol
end
