#=
This file contains all the helper functions called by the main functions to run Castawave. This includes functions for mathematical computations, such as the conformal mapping or timestepping scheme, and functions for additional tasks such as output writing and visualization.
=#

using DrWatson
@quickactivate "Castawave"
using Plots
using LinearAlgebra 

include("Constants.jl")

function conformalMap(R::Vector)
    #=
    conformalMap is a function that takes complex values R(ξ) and conformally transforms them. It is assumed that ξ is the Lagrangian complex spatial coordinate, where R is the complex surface.

    See section 3 in Longuet-Higgins and Cokelet (1976) or Dold (1992) equation 4.6 for details.

    Input:
    R - A complex vector representing the X + iY positions of each particle on the surface 
    period - The periodicity of the wave in physical units

    Output:
    Ω - conformally mapped closed contour of the fluid surface
    r - the modulus of Ω
    θ - the phase angle of Ω
    =#

    Ω = exp.(- im * R);

    r = abs.(Ω)
    θ = angle.(Ω)

    return Ω, r, θ
end

function conformalDepth(h)
    #=
    conformalDepth is a function that takes the real depth h and transforms it to the more useful conformal value H, which tends to 0 at inifinite depth. To avoid working with infinite values, the infinite depth assumption is taken when the user inputs 0 as the real h.
    =#
    if h == 0.0
        H = 0.0
    else
        H = exp(-2 * h)
    end

    return H 
end

#=
The following DDI1 and DDI2 functions are used for taking the first and second order tangential derivatives with respect to the particle label ξ. They do this according to Dold's weighted coefficients method for Lagrangian polynomial interpolation. They are implemented in the main function with respect to conformal methods such that no offset or halo boundary is needed to deal with periodic boundary conditions.ive.
=#
function DDI1(Ω::Vector, N, q=0)
    #=
    DDI1 is a function that implements the 11-point Lagrangian polynomial interpolation for the first coefficient from eq... , therefore returning the first order derivative value at the center point.
    =#
    if q == 0
        Ω_p = zeros(Complex, N)
    elseif q == 1
        Ω_p = zeros(N)
    end

    for i in 1:N
        points = [Ω[mod1(i+1, N)] - Ω[mod1(i-1, N)], Ω[mod1(i+2, N)] - Ω[mod1(i-2, N)], Ω[mod1(i+3, N)] - Ω[mod1(i-3, N)], Ω[mod1(i+4, N)] - Ω[mod1(i-4, N)], Ω[mod1(i+5, N)] - Ω[mod1(i-5, N)]]
        Ω_p[i] = sum(COEFFICIENTS1 .* points)
    end

    return Ω_p
end

function DDI2(Ω::Vector, N, q=0)
    if q == 0
        Ω_pp = zeros(Complex, N)
    elseif q == 1
        Ω_pp = zeros(N)
    end

    for i in 1:N
        points = [Ω[i], Ω[mod1(i+1, N)] + Ω[mod1(i-1, N)], Ω[mod1(i+2, N)] + Ω[mod1(i-2, N)], Ω[mod1(i+3, N)] + Ω[mod1(i-3, N)], Ω[mod1(i+4, N)] + Ω[mod1(i-4, N)], Ω[mod1(i+5, N)] + Ω[mod1(i-5, N)]]
        Ω_pp[i] = sum(COEFFICIENTS2 .* points)
    end

    return Ω_pp
end

function ABMatrices(Ω, Ω_ξ, Ω_ξξ, N, H=0.0)
    #=
    The ABMatrices function sets up the A and B matrices from Dold eq. 4.13, using a necesary conditional double for loop. While it is a computationally expensive function, it is separated from the matrix inversion step because it is only needed once per timestep.

    Input:
    Ω - complex vector of conformally mapped positions of particles
    Ω_ξ - first-order derivative of Ω with respect to particle labels ξ
    N - number of particles
    H - parameter specified in Dold for bottom depth boundary condition. Defaulted to 0 for infinite depth.

    Output:
    A - matrix
    B - matrix
    ℵ - matrix of π - B for use in NormalInversion
    =#

    C = zeros(Complex, N, N)    # C = A + iB

    # Does this have to be a double loop because of condition? Access column-first for optimization. The column index j corresponds to ξ prime, the rows to regular ξ.
    # How to optimize the H=0 condition? Has to be in the loops to avoid another long double loop, but need it be checked each time?

    if H == 0.0
        for ξ_p in 1:N
            for ξ in 1:N
                if ξ_p == ξ     
                    C[ξ,ξ_p] = Ω_ξξ[ξ] / (2 * Ω_ξ[ξ])
                else
                    C[ξ,ξ_p] = Ω_ξ[ξ] / (Ω[ξ] - Ω[ξ_p])
                end
            end
        end
    else
        for ξ_p in 1:N
            for ξ in 1:N
                C[ξ,ξ_p] = -conj((H * Ω_ξ[ξ] / conj(Ω[ξ_p])) / (Ω[ξ] * (Ω[ξ] - H / conj(Ω[ξ_p]))))
                if ξ_p == ξ     
                    C[ξ,ξ_p] += Ω_ξξ[ξ] / (2 * Ω_ξ[ξ])
                else
                    C[ξ,ξ_p] += Ω_ξ[ξ] / (Ω[ξ] - Ω[ξ_p])
                end
            end
        end
    end
        
    A = real(C) 
    B = imag(C)

    ℵ = π * I - B   # Identity matrix I from LinearAlgebra to subtract B from pi diagonal.

    return A,B,ℵ
end

function NormalInversion(ϕ, A, ℵ, N)
    #= 
    NormalInversion is a function that implements the matrix inversion method from Dold eq 4.13 in order to compute the normal derivative of the scalar velocity potential. The method is based on using a conformal mapping and the Cauchy integral theorem for solving the Laplacian equation. The subtlety lies in the issue that, for solely surface particles at b=0, there is no Lagrangian normal derivative as there are no particles above or below. For efficiency, due to the need of the tangential derivative, it is first computed and returned along with the normal derivative here.
        
    The implementation here is mathematically simplified from Dold's formula to take the form of an Ax = b matrix equation, and is thus quite compact and optimized.

    Input:
    ϕ - real vector of scalar velocity potential
    A - 
    ℵ - 
    N - number of particles

    Output:
    ϕ_ξ - real vector of tangential partial derivative of with respect to particle label ξ
    ϕ_ν - real vector of normal partial derivative scaled by
    =#

    ϕ_ξ = DDI1(ϕ, N, 1)
    ϕ_ξξ = DDI2(ϕ, N, 1)

    # Important: here the * is not element wise to get the sum A*ϕ_ξ for each one-element row entry of the resulting column vector, while the difference is element wise to subtract ϕ_ξξ[i] from each of the summed entries.
    b = (A * ϕ_ξ) .- ϕ_ξξ    

    # Ax = b using the efficient \ operator, where x is the vector of tangential derivatives
    ϕ_ν = ℵ \ b

    return ϕ_ξ, ϕ_ν
end

function PhiTimeDer(R_ξ, ϕ_ξ, ϕ_ν, Y, L)
    #=
    PhiTimeDer is a function that calculates the Lagrangian and Eulerian time derivative of ϕ from the dynamic Bernoulli condition. 
        
    Dold's method of timestepping involves a Taylor series with multiple order of time derivatives of X, Y, ϕ, which are found in Appendix A from differentiating the Bernoulli equation. This requires a subtle procedure of recomputing the tangential and normal derivative of a higher time derivative of ϕ from the previous order. 
    For alternative (valid) timestepping schemes, notably RK4, only the first derivative, found directly from the Bernoulli equation, is needed. 

    Input:
    R_ξ - 
    ϕ_ξ - 
    ϕ_ν - 
    Y - 
    
    Output:
    ϕ_D - Material Lagrangian time derivative of the velocity potential
    ϕ_t - Partial Eulerian time derivative of the velocity potential
    =#
    ϕ_D = 0.5 .* (ϕ_ξ.^2 .+ ϕ_ν.^2) ./ abs.(R_ξ).^2 .- GRAVITY .* Y
    ϕ_t = -0.5 .* (ϕ_ξ.^2 .+ ϕ_ν.^2) ./ abs.(R_ξ).^2 .- GRAVITY .* Y

    return ϕ_D, ϕ_t
end



function RK4i(dt, f::Function, N, X, Y, ϕ, L, H,smoothing)
    #=
    RK4 is a function that implements the Runge-Kutta fourth order timestepping scheme for a triple of vectors (X, Y, ϕ)

    Input:
    dt - timestep (fixed for now)
    f - evaluation function of the dy/dx form
    N - number of particles
    X - real vector of particle x-positions on the surface
    Y - real vector of particle y-positions on the surface
    ϕ - real vector of scalar velocity potential for particles on the surface
    L - Length of domain (meters )
    H - Depth of domain (defaults to zero)
    =#

    k1X, k1Y, k1ϕ = f(N, X, Y, ϕ, L, H,smoothing)
    k2X, k2Y, k2ϕ = f(N, X .+ dt ./ 2 .* k1X, Y .+ dt ./ 2 .* k1Y, ϕ .+ dt ./ 2 .* k1ϕ, L, H,smoothing)
    k3X, k3Y, k3ϕ = f(N, X .+ dt ./ 2 .* k2X, Y .+ dt ./ 2 .* k2Y, ϕ .+ dt ./ 2 .* k2ϕ, L, H,smoothing)
    k4X, k4Y, k4ϕ = f(N, X .+ dt .* k3X, Y .+ dt .* k3Y, ϕ .+ dt .* k3ϕ, L, H,smoothing)
    
    Xn = X .+ dt ./ 6 .* (k1X .+ 2 .* k2X .+ 2 .* k3X .+ k4X)     
    Yn = Y .+ dt ./ 6 .* (k1Y .+ 2 .* k2Y .+ 2 .* k3Y .+ k4Y)     
    ϕn = ϕ .+ dt ./ 6 .* (k1ϕ .+ 2 .* k2ϕ .+ 2 .* k3ϕ .+ k4ϕ)

    return Xn, Yn, ϕn
end

function RealPhi(R_ξ, ϕ_ξ, ϕ_ν)
    #=
    RealPhi is a small helper function that is used to transform back from the conformally mapped variables to the real plane when timestepping the real X, Y, ϕ arrays. Specically, the relation formula between the complex ϕ_x + iϕ_y and the ξ and ν partial derivatives is exploited.

    Input:
    R_ξ - complex vector representation of position vectors X + iY
    ϕ_ξ - real vector of partial derivatives of ϕ with respect to label ξ
    ϕ_ν - real vector of normal partial derivatives of ϕ scaled by 

    Output:
    ϕ_x - real vector of U velocity for particles on the surface
    ϕ_y - real vector of V velocity for particles on the surface
    =#

    ϕ_grad = ((ϕ_ξ) .+ im .* (ϕ_ν)) ./ (conj.(R_ξ))
    ϕ_x = real.(ϕ_grad)
    ϕ_y = imag.(ϕ_grad)

    return ϕ_x, ϕ_y
end

function TaylorTimestep(dt, f, f_t = 0, f_tt = 0, f_ttt = 0)
    # We give a function of ζ and t, ϕ or r (x, y). Within a timestep, time is fixed so the functions given are just spatial. Then g = f(t+dt) by truncated Taylor expansion.
    g = f .+ dt .* f_t .+ 0.5 .* dt^2 .* f_tt .+ (1/6) .* dt^3 .* f_ttt
    return g
end

function TaylorTimestepBD(dt, f, f_t, f_tt, fbd1, fbd2, fbd3)
    # We give a function of ζ and t, ϕ or r (x, y). Within a timestep, time is fixed so the functions given are just spatial. Then g = f(t+dt) by truncated Taylor expansion.
    g = f .+ dt .* f_t .+ 0.5 .* dt^2 .* f_tt .+ (1/24) .* dt^4 .* fbd1 .+ (1/60) .* dt^5 .* fbd2 .+ (1/120) .* dt^5 .* fbd3
    return g
end

function mwl(X_ξ, Y)
    # mwl is a simple function that computes the mean water level.
    return 1 / (2 * pi) * sum(Y .* X_ξ)
end

function smooth(N, Ω, q=1)
    println("smooth")
    if q == 0
        Ω_sm = zeros(Complex, N)
    elseif q == 1
        Ω_sm = zeros(N)
    end

    for i in 1:N

        if i == 4 && false
            points = [Ω[i], Ω[mod1(i+1, N)] + Ω[mod1(i-1, N)], Ω[mod1(i+2, N)] + Ω[mod1(i-2, N)], Ω[mod1(i+3, N)] + Ω[mod1(i-3, N)], Ω[mod1(i+4, N)] + Ω[mod1(i-4, N)], Ω[mod1(i+5, N)] + Ω[mod1(i-5, N)] + 2*π, Ω[mod1(i+6, N)] + Ω[mod1(i-6, N)], Ω[mod1(i+7, N)] + Ω[mod1(i-7, N)]]
        else
            points = [Ω[i], Ω[mod1(i+1, N)] + Ω[mod1(i-1, N)], Ω[mod1(i+2, N)] + Ω[mod1(i-2, N)], Ω[mod1(i+3, N)] + Ω[mod1(i-3, N)], Ω[mod1(i+4, N)] + Ω[mod1(i-4, N)], Ω[mod1(i+5, N)] + Ω[mod1(i-5, N)], Ω[mod1(i+6, N)] + Ω[mod1(i-6, N)], Ω[mod1(i+7, N)] + Ω[mod1(i-7, N)]]
        end

        δM = sum(SMOOTHCOEFFICIENTS .* points) / (2^14)

        #Ω_sm[i] = Ω[i] - δM
        if iseven(i)
            Ω_sm[i] = Ω[i] - δM
        else
            Ω_sm[i] = Ω[i] + δM
        end
    end

    return Ω_sm
end

function simpsons_rule_periodic(X, Y)
    n = length(X)
    h = diff(X)
    sum = 0   
    for i in 1:2:n-3
        h1 = h[i]
        h2 = h[i+1]
        hph = h2 + h1
        hdh = h2/h1
        hmh = h2*h1
        sum += (hph/6)*((2 - hdh)*Y[i] + (hph^2/hmh)*Y[i+1] + (2 - 1/hdh)*Y[i+2])
    end
    # handle last point periodically
    i = n-1
    h1 = h[i]
    h2 = h[1]
    hph = h2 + h1
    hdh = h2/h1
    hmh = h2*h1
    sum += (hph/6)*((2 - hdh)*Y[i] + (hph^2/hmh)*Y[i+1] + (2 - 1/hdh)*Y[1])
    
    return sum
end


function computeEnergy(sol,n,Δt,tf)
    energy = []
    MWL_check = []
    for t ∈ 0:Δt:tf
        x = sol(t)[1:n]
        y = sol(t)[n+1:2*n]
        ẋ = sol(t,Val{1})[1:n]
        ẏ = sol(t,Val{1})[n+1:2*n]
        ϕ = sol(t)[2*n+1:end]
        R = x .+ im.*y
        Ω,r,θ = conformalMap(R)
        Ω_ξ = DDI1(Ω, n)
        R_ξ = (im ./ Ω) .* Ω_ξ
        q = ẋ .+ im.*ẏ
        ϕ_n = imag.(q .* abs.(R_ξ) ./ R_ξ)
        # Other way of getting KE from Balk 
        B_ξ = ẋ.*imag.(R_ξ) .- ẏ.*real.(R_ξ)
        integrand = -1/2 * ϕ .* B_ξ
        #KE = simpsons_rule_periodic(x,ϕ.* ϕ_n/2)
        KE = simpsons_rule_periodic(1:n,integrand)
        #println(KE)
        PE = simpsons_rule_periodic(1:n,GRAVITY/2 * (y).^2 .* real.(R_ξ))
        append!(MWL_check,simpsons_rule_periodic(x,y)/(2π)) # Eulerian MWL
        #append!(MWL_check,sum(y)/n) # Lagrangian MWL
        append!(energy,KE + PE)
    end
    return energy, MWL_check
end
