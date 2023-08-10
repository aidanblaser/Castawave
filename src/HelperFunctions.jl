#=
This file contains all the helper functions called by the main functions to run Castawave. This includes functions for mathematical computations, such as the conformal mapping or timestepping scheme, and functions for additional tasks such as output writing and visualization.
=#

using DrWatson
@quickactivate "PacketDrift"
using Plots
using LinearAlgebra 

function conformalMap(R::Vector)
    #=
    conformalMap is a function that takes complex values R(ξ) and conformally transforms them. It is assumed that ξ is the Lagrangian complex spatial coordinate, where R is the complex surface.

    See section 3 in Longuet-Higgins and Cokelet (1976) or Dold (1992) equation 4.6 for details.

    Input:
    R - A complex vector representing the X + iY positions of each particle on the surface 

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

#=
The following functions take the Lagrangian derivative with the condition of the periodicity implemented as an optional parameter. When chosen, the derivative is taken with respect to the X - a offset values, thus taking the periodic difference function's derivative.
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

#=
function offset(Ω, N)
    ω = zeros(Complex, N)
    for j in 1:N
        ω[j] = Ω[j] - (j * 2 * π / N)
    end
    return ω
end

function DDI1full(Ω::Vector, N, p=0, q=0)
    #=
    TenthOrderFirstD is a function that finds the first-order derivative for all points of a given vector. The derivatives are calculated at the center point of an 11-point Lagrangian interpolation ploynomial. The method and coefficients used are described in Dold Appendix C. The parameters assume the inputted vector is periodic such that the 11 interpolation points can be taken as if from a circular array.

    Input: 
    Ω - Array of points Ω 
    N - Number of points
    
    Output:
    Ω_p - Vector of length N of first derivatives.
    =#
    if q == 0
        Ω_p = zeros(Complex, N)
    elseif q == 1
        Ω_p = zeros(N)
    end
    if p == 0
        for i in 1:N
            pointsm = [Ω[mod1(i+1, N)] - Ω[mod1(i-1, N)], Ω[mod1(i+2, N)] - Ω[mod1(i-2, N)], Ω[mod1(i+3, N)] - Ω[mod1(i-3, N)], Ω[mod1(i+4, N)] - Ω[mod1(i-4, N)], Ω[mod1(i+5, N)] - Ω[mod1(i-5, N)]]
            Ω_p[i] = (WEIGHTS1 \ (FACTOR1 * pointsm))[1]
        end
    elseif p == 1
        ω = offset(Ω, N) 
        for i in 1:N
            pointsm = [ω[mod1(i+1, N)] - ω[mod1(i-1, N)], ω[mod1(i+2, N)] - ω[mod1(i-2, N)], ω[mod1(i+3, N)] - ω[mod1(i-3, N)], ω[mod1(i+4, N)] - ω[mod1(i-4, N)], ω[mod1(i+5, N)] - ω[mod1(i-5, N)]]
            Ω_p[i] = (WEIGHTS1 \ (FACTOR1 * pointsm))[1]
        end
    elseif p == 2
        # l = real(Ω[N]) - real(Ω[1]) + (real(Ω[2]) - real(Ω[1])) 
        l = 2*π
        q = zeros(Complex, N+10)
        for i in 1:5
            q[i] = real(Ω[mod1(i-5, N)]) - l + imag(Ω[mod1(i-5, N)])im
        end
        for i in 6:N+5
            q[i] = Ω[i-5]
        end
        for i in N+6:N+10
            q[i] = real(Ω[mod1(i-5, N)]) + l + imag(Ω[mod1(i-5, N)])im
        end
        for i in 6:N+5
            pointsm = [q[i+1] - q[i-1], q[i+2] - q[i-2], q[i+3] - q[i-3], q[i+4] - q[i-4], q[i+5] - q[i-5]]
            Ω_p[i-5] = (WEIGHTS1 \ (FACTOR1 * pointsm))[1]
        end
    elseif p == 3
        # l = real(Ω[N]) - real(Ω[1]) + (real(Ω[2]) - real(Ω[1])) 
        l = 2*π
        q = zeros(Complex, N+10)
        for i in 1:5
            q[i] = Ω[mod1(i-5, N)] - l
        end
        for i in 6:N+5
            q[i] = Ω[i-5]
        end
        for i in N+6:N+10
            q[i] = Ω[mod1(i-5, N)] + l 
        end
        for i in 6:N+5
            pointsm = [q[i+1] - q[i-1], q[i+2] - q[i-2], q[i+3] - q[i-3], q[i+4] - q[i-4], q[i+5] - q[i-5]]
            Ω_p[i-5] = (WEIGHTS1 \ (FACTOR1 * pointsm))[1]
        end
    end

    return Ω_p
end

function DDI2full(Ω::Vector, N, p=0, q=0)
    #=
    TenthOrderFirstD is a function that finds the second-order derivative for all points of a given vector. The derivatives are calculated at the center point of an 11-point Lagrangian interpolation ploynomial. The method and coefficients used are described in Dold Appendix C. The parameters assume the inputted vector is periodic such that the 11 interpolation points can be taken as if from a circular array.

    Input: 
    Ω - Array of points Ω 
    N - Number of points
    
    Output:
    Ω_p - Vector of length N of second derivatives.
    =#

    if q == 0
        Ω_pp = zeros(Complex, N)
    elseif q == 1
        Ω_pp = zeros(N)
    end
        
    if p == 0
        for i in 1:N
            pointsm = [Ω[i], Ω[mod1(i+1, N)] + Ω[mod1(i-1, N)], Ω[mod1(i+2, N)] + Ω[mod1(i-2, N)], Ω[mod1(i+3, N)] + Ω[mod1(i-3, N)], Ω[mod1(i+4, N)] + Ω[mod1(i-4, N)], Ω[mod1(i+5, N)] + Ω[mod1(i-5, N)]]
            Ω_pp[i] = (WEIGHTS2 \ (FACTOR2 * pointsm))[1]
        end
    elseif p == 1
        ω = offset(Ω, N)
        for i in 1:N
            pointsm = [ω[i], ω[mod1(i+1, N)] + ω[mod1(i-1, N)], ω[mod1(i+2, N)] + ω[mod1(i-2, N)], ω[mod1(i+3, N)] + ω[mod1(i-3, N)], ω[mod1(i+4, N)] + ω[mod1(i-4, N)], ω[mod1(i+5, N)] + ω[mod1(i-5, N)]]
            Ω_pp[i] = (WEIGHTS2 \ (FACTOR2 * pointsm))[1]
        end
    end
    
    return 2 * Ω_pp
end
=#

#=
#LAGRANGIAN 11 POINT INTERPOLATION FOR TANGENTIAL DERIVATIVES (C.5 formula). The input is a vector for the N particles specifying a property field. The output is the complex derivative Vector.
function TenthOrderFirstD(Ω::Vector, m=0, N=100)
    # Initialize the array of derivative values.
    Ω_p = zeros(Complex, N)

    if m == 0
        for i in 1:N
            pointsm = [Ω[mod1(i+1, N)] - Ω[mod1(i-1, N)], Ω[mod1(i+2, N)] - Ω[mod1(i-2, N)], Ω[mod1(i+3, N)] - Ω[mod1(i-3, N)], Ω[mod1(i+4, N)] - Ω[mod1(i-4, N)], Ω[mod1(i+5, N)] - Ω[mod1(i-5, N)]]
            Ω_p[i] = (WEIGHTS1 \ (FACTOR1 * pointsm))[1]
        end
    elseif m == 1
        # l = real(Ω[N]) - real(Ω[1]) + (real(Ω[2]) - real(Ω[1])) 
        l = 1.0
        q = zeros(Complex, N+10)
        for i in 1:5
            q[i] = real(Ω[mod1(i-5, N)]) - l + imag(Ω[mod1(i-5, N)])im
        end
        for i in 6:N+5
            q[i] = Ω[i-5]
        end
        for i in N+6:N+10
            q[i] = real(Ω[mod1(i-5, N)]) + l + imag(Ω[mod1(i-5, N)])im
        end
        for i in 6:N+5
            pointsm = [q[i+1] - q[i-1], q[i+2] - q[i-2], q[i+3] - q[i-3], q[i+4] - q[i-4], q[i+5] - q[i-5]]
            Ω_p[i-5] = (WEIGHTS1 \ (FACTOR1 * pointsm))[1]
        end
    elseif m == 2
        l = 1.0
        for i in 1:5
            pointsm = [Ω[i+1] - modval(Ω[mod1(i-1, N)], l, -1), Ω[i+2] - modval(Ω[mod1(i-2, N)], l, -1), Ω[i+3] - modval(Ω[mod1(i-3, N)], l, -1), Ω[i+4] - modval(Ω[mod1(i-4, N)], l, -1), Ω[i+5] - modval(Ω[mod1(i-5, N)], l, -1)]
            Ω_p[i] = (WEIGHTS1 \ (FACTOR1 * pointsm))[1]
        end
        for i in 6:N-5
            pointsm = [Ω[i+1] - Ω[i-1], Ω[i+2] - Ω[i-2], Ω[i+3] - Ω[i-3], Ω[i+4] - Ω[i-4], Ω[i+5] - Ω[i-5]]
            Ω_p[i] = (WEIGHTS1 \ (FACTOR1 * pointsm))[1]
        end
        for i in N-4:N
            pointsm = [modval(Ω[mod1(i+1, N)], l, 1) - Ω[i-1], modval(Ω[mod1(i+2, N)], l, 1) - Ω[i-2], modval(Ω[mod1(i+3, N)], l, 1) - Ω[i-3], modval(Ω[mod1(i+4, N)], l, 1) - Ω[i-4], modval(Ω[mod1(i+5, N)], l, 1) - Ω[i-5]]
            Ω_p[i] = (WEIGHTS1 \ (FACTOR1 * pointsm))[1]
        end
    end
    return Ω_p
end

function modval(mods, l, p)
    if p == 1
        md = real(mods) + l + imag(mods)im
    elseif p == -1
        md = real(mods) - l + imag(mods)im
    end
    return md     
end

function TenthOrderSecondD(Ω::Vector)
    Ω_pp = zeros(Complex, N)

    for i in 1:N
        pointsm = [Ω[i], Ω[mod1(i+1, N)] + Ω[mod1(i-1, N)], Ω[mod1(i+2, N)] + Ω[mod1(i-2, N)], Ω[mod1(i+3, N)] + Ω[mod1(i-3, N)], Ω[mod1(i+4, N)] + Ω[mod1(i-4, N)], Ω[mod1(i+5, N)] + Ω[mod1(i-5, N)]]
        Ω_pp[i] = (WEIGHTS2 \ (FACTOR2 * pointsm))[1]
    end

    return 2 * Ω_pp
end

function TenthOrderFirstDIndividual(nd::Vector)
    #In this function, note the 11 point nodes from 1:11 are actually from -5 to 5 in the paper, where the sixth index is the midpoint.
    factor1 = 
    [2100 -600 150 -25 2; 
    -70098 52428 -14607 2522 -205;
    1938 -1872 783 -152 13;
    -378 408 -207 52 -5;
    42 -48 27 -8 1]

    w1 = [2520, 181440, 34560, 120960, 725760]
    m1 = inv(diagm(w1))
    c1 = [nd[7] - nd[5], nd[8] - nd[4], nd[9] - nd[3], nd[10] - nd[2], nd[11] - nd[1]]

    coeffs = m1 * factor1 * c1  #this gives the odd coefficients
    return coeffs[1]
end

function TenthOrderSecondDIndividual(nd::Vector)
    factor2 = 
    [-73766 42000 -6000 1000 -125 8;
    192654 -140196 52428 -9738 1261 -82;
    -12276 9690 -4680 1305 -190 13;
    462 -378 204 -69 13 -1;
    -252 210 -120 45 -10 1]

    w2 = [50400, 362880, 172800, 120960, 3628800]
    m2 = inv(diagm(w2))
    c2 = [nd[6], nd[7] + nd[5], nd[8] + nd[4], nd[9] + nd[3], nd[10] + nd[2], nd[11] + nd[1]]

    coeffs = m2 * factor2 * c2  #this gives the even coefficients
    return 2 * coeffs[1]
end
=#

function ABMatrices(Ω, Ω_ξ, Ω_ξξ, N, H=0)
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

    if H == 0
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
                C[i,j] = -conj((H * Ω_ξ[ξ] / conj(Ω[ξ_p])) / (Ω[ξ] * (Ω[ξ] - H / conj(Ω[ξ_p]))))
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

function PhiTimeDer(R_ξ, ϕ_ξ, ϕ_ν, Y)
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

    ϕ_D = 0.5 .* (ϕ_ξ.^2 .+ ϕ_ν.^2) ./ abs.(R_ξ).^2 .- 9.81 .* Y
    ϕ_t = -0.5 .* (ϕ_ξ.^2 .+ ϕ_ν.^2) ./ abs.(R_ξ).^2 .- 9.81 .* Y

    return ϕ_D, ϕ_t
end

#=
function PhiSecondTimeDer(R_ξ, ϕ_x, ϕ_y, ϕ_t, A, ℵ, N, ϕ_ξ, ϕ_ν)
    X_ξ = real.(R_ξ)
    Y_ξ = imag.(R_ξ)

    U_ξ = DDI1(ϕ_x, N)
    V_ξ = DDI1(ϕ_y, N)

    ϕ_tξ, ϕ_tν = NormalInversion(ϕ_t, A, ℵ, N)
    ∇ϕ_t = (ϕ_tξ .+ im .* ϕ_tν) ./ conj.(R_ξ)
    # ∇ϕ_t = (ϕ_tξ .+ im .* ϕ_tν) ./ conj.(R_ξ) .- (U_ξ .- im .* V_ξ) .* (ϕ_ξ .+ im .* ϕ_ν) ./ conj.(R_ξ).^2
    ϕ_xt = real.(∇ϕ_t)
    ϕ_yt = imag.(∇ϕ_t)

    ux = (U_ξ .* X_ξ .- V_ξ .* Y_ξ) ./ (abs.(R_ξ).^2)
    vy = -ux
    vx = (U_ξ .* Y_ξ .+ V_ξ .* X_ξ) ./ (abs.(R_ξ).^2)
    uy = vx

    U_D = ϕ_xt .+ ϕ_x .* ux .+ ϕ_y .* uy
    V_D = ϕ_yt .+ ϕ_x .* vx .+ ϕ_y .* vy
    # DUDT = ϕ_xt .+ (ϕ_x .* ux .+ ϕ_y .* vy) .* ϕ_x
    # DVDT = ϕ_yt .+ (ϕ_x .* vx .+ ϕ_y .* vy) .* ϕ_y
    
    ϕ_DD = -9.81 .* ϕ_y .+ ϕ_x .* U_D .+ ϕ_y .* V_D

    return U_D, V_D, ϕ_DD
end

function PhiThirdTimeDer(R_ξ, ϕ_x, ϕ_y, ϕ_t, A, ℵ, N, ϕ_ξ, ϕ_ν)
    X_ξ = real.(R_ξ)
    Y_ξ = imag.(R_ξ)
    # U_ξ = DDI1(real.((ϕ_ξ .+ im .* ϕ_ν) ./ conj.(R_ξ)), N)
    # V_ξ = DDI1(imag.((ϕ_ξ .+ im .* ϕ_ν) ./ conj.(R_ξ)), N)
    U_ξ = DDI1(ϕ_x, N)
    V_ξ = DDI1(ϕ_y, N)

    ϕ_tξ, ϕ_tν = NormalInversion(ϕ_t, A, ℵ, N)

    ∇ϕ_t = (ϕ_tξ .+ im .* ϕ_tν) ./ conj.(R_ξ) .- (U_ξ .- im .* V_ξ) .* (ϕ_ξ .+ im .* ϕ_ν) ./ conj.(R_ξ).^2
    ϕ_xt = real.(∇ϕ_t)
    ϕ_yt = imag.(∇ϕ_t)

    ux = (U_ξ .* X_ξ .- V_ξ .* Y_ξ) ./ (abs.(R_ξ).^2)
    vy = -ux
    vx = (U_ξ .* Y_ξ .+ V_ξ .* X_ξ) ./ (abs.(R_ξ).^2)
    uy = vx

    U_D = ϕ_xt .+ ϕ_x .* ux .+ ϕ_y .* uy
    V_D = ϕ_yt .+ ϕ_x .* vx .+ ϕ_y .* vy
    
    ϕ_DD = -9.81 .* ϕ_y .+ ϕ_x .* U_D .+ ϕ_y .* V_D

    U_ξt = DDI1(ϕ_xt, N)
    V_ξt = DDI1(ϕ_yt, N)

    ϕ_tt = - (9.81 .* ϕ_y .+ ϕ_x .* U_D .+ ϕ_y .* V_D .+ ϕ_x .* ϕ_xt .+ ϕ_y .* ϕ_yt)
    
    ϕ_ttξ, ϕ_ttν = NormalInversion(ϕ_tt, A, ℵ, N)

    ∇ϕ_tt = (ϕ_ttξ .+ im .* ϕ_ttν) ./ conj.(R_ξ) .- 2 .* (U_ξ .- im .* V_ξ) .* (ϕ_tξ .+ im .* ϕ_tν) ./ conj.(R_ξ).^2 .- (U_ξt .- im .* V_ξt) .* (ϕ_ξ .+ im .* ϕ_ν) ./ conj.(R_ξ).^2 .+ 2 .* ((U_ξ .- im .* V_ξ) .^ 2) .* (ϕ_ξ .+ im .* ϕ_ν) ./ conj.(R_ξ).^3

    ϕ_xtt = real.(∇ϕ_tt)
    ϕ_ytt = imag.(∇ϕ_tt)

    uxt = 1
    vyt = 1
    vxt = 1
    uyt = 1

    uxx = 1
    vxx = 1

    U_DD = ϕ_xtt .+ 2 .* (ϕ_x .* uxt .+ ϕ_y .* uyt) .+ (ϕ_xt .* ux .+ ϕ_yt .* vx) .+ ((ux.^2 .+ vx.^2) .* ϕ_x .+ (ϕ_x.^2 - ϕ_y.^2) .* uxx .+ (2 .* ϕ_x .* ϕ_y .* vxx))
    V_DD = ϕ_ytt .+ 2 .* (ϕ_x .* vxt .+ ϕ_y .* vyt) .+ (ϕ_xt .* vx .- ϕ_yt .* ux) .+ ((ux.^2 .+ vx.^2) .* ϕ_y .+ (ϕ_x.^2 - ϕ_y.^2) .* vxx .- (2 .* ϕ_x .* ϕ_y .* uxx))



    return U_D, V_D, ϕ_DD
end
=#


function RK4i(dt, f::Function, N, X, Y, ϕ)
    #=
    RK4 is a function that implements the Runge-Kutta fourth order timestepping scheme for a triple of vectors (X, Y, ϕ)

    Input:
    dt - timestep (fixed for now)
    f - evaluation function of the dy/dx form
    N - number of particles
    X - real vector of particle x-positions on the surface
    Y - real vector of particle y-positions on the surface
    ϕ - real vector of scalar velocity potential for particles on the surface
    =#

    k1X, k1Y, k1ϕ = f(N, X, Y, ϕ)
    k2X, k2Y, k2ϕ = f(N, X .+ dt ./ 2 .* k1X, Y .+ dt ./ 2 .* k1Y, ϕ .+ dt ./ 2 .* k1ϕ)
    k3X, k3Y, k3ϕ = f(N, X .+ dt ./ 2 .* k2X, Y .+ dt ./ 2 .* k2Y, ϕ .+ dt ./ 2 .* k2ϕ)
    k4X, k4Y, k4ϕ = f(N, X .+ dt .* k3X, Y .+ dt .* k3Y, ϕ .+ dt .* k3ϕ)
    
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

function smooth(N, Ω, q=0)
    if q == 0
        Ωsm = zeros(Complex, N)
    elseif q == 1
        Ωsm = zeros(N)
    end

    for i in 1:N
        points = [Ω[i], Ω[mod1(i+1, N)] + Ω[mod1(i-1, N)], Ω[mod1(i+2, N)] + Ω[mod1(i-2, N)], Ω[mod1(i+3, N)] + Ω[mod1(i-3, N)], Ω[mod1(i+4, N)] + Ω[mod1(i-4, N)], Ω[mod1(i+5, N)] + Ω[mod1(i-5, N)], Ω[mod1(i+6, N)] + Ω[mod1(i-6, N)], Ω[mod1(i+7, N)] + Ω[mod1(i-7, N)]]

        δM = sum(SMOOTHCOEFFICIENTS .* points) / (2^14)

        Ω_sm[i] = Ω[i] - δM
        # if iseven(i)
        #     Ω_sm[i] = Ω[i] - δM
        # else
        #     Ω_sm[i] = Ω[i] + δM
        # end
    end

    return Ωsm
end

function turningAngle(N, Ω)
    ta = zeros(N)

    for i in 1:N
        ta[i] = abs(angle(Ω[mod1(i+1, N)] - Ω[i]) - angle(Ω[i] - Ω[mod1(i-1, N)]))
    end

    return maximum(ta)
end

function rms()
    1
end


#=
#EXTRA FUNCTIONS FOR NOW
# Transforming conformal tangential derivatives to real plane, from Dold code.
function GradientSTransform(u, v, ϕ_u, ϕ_v)
    ϕ_x = (ϕ_u .* v .- ϕ_v .* u) ./ (ϕ_u .* ϕ_u .+ ϕ_v .* ϕ_v)
    ϕ_y = (ϕ_v .* v .+ ϕ_u .* u) ./ (ϕ_u .* ϕ_u .+ ϕ_v .* ϕ_v)
    return ϕ_x, ϕ_y
end

# For potential future use, implements the cosine transform for Chebyshev spacing
function ChebyshevMap(n::Int64)
    l = (n - 1) / 2
    c = zeros(n)
    for i in 1:n
        g = π / 2
        f = 2 * i + 1
        h = (n + 1)
        c[i] = - cos(g * f / h)
    end
end

# Simple function that creates the gradient values as a Vector of Tuples of the partial derivatives.
function GradientST(ϕ_ξ, ϕ_ν, N)
    ∇ϕ = [(ϕ_ξ[i], ϕ_ν[i]) for i in 1:N]
    return ∇ϕ
end

# function Phi2Time(R_ξ, ϕ_x, ϕ_y, ϕ_ξ, ϕ_ν, ϕ_t, A, ℵ, N)
#     X_ξ = real.(R_ξ)
#     Y_ξ = imag.(R_ξ)
#     U_ξ = DDI1(ϕ_x, N)
#     V_ξ = DDI1(ϕ_y, N)

#     ϕ_tξ, ϕ_tν = NormalInversion(ϕ_t, A, ℵ, N)
#     ϕ_tt = (ϕ_ξ .* ϕ_tξ .+ ϕ_ν .* ϕ_tν) ./ (abs.(R_ξ).^2) .- (ϕ_ξ .^ 2 .* (X_ξ .* U_ξ .+ Y_ξ .* V_ξ) .+ ϕ_ν .^ 2 .* (X_ξ .* U_ξ .+ Y_ξ .* V_ξ)) ./ (abs.(R_ξ).^4) - 9.81 .* ϕ_y
    
#     return ϕ_tt
# end

=#