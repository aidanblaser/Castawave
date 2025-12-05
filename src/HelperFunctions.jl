#=
This file contains all the helper functions called by the main functions to run Castawave. This includes functions for mathematical computations, such as the conformal mapping or timestepping scheme, and functions for additional tasks such as output writing and visualization.
=#

using DrWatson
@quickactivate "Castawave"
using LinearAlgebra 
using Statistics
using Polynomials

include("Constants.jl")

#= TODO
- Properly broadcast creation of A,B matrices for finite depth 
- Try spectral derivatives in DDI1 and DDI2
=#

function conformalMap(R::AbstractVector{<:Number})
    #=
    conformalMap is a function that takes complex values R(ξ) and conformally transforms them. It is assumed that ξ is the Lagrangian complex spatial coordinate,
    where R is the complex surface.

    See section 3 in Longuet-Higgins and Cokelet (1976) or Dold (1992) equation 4.6 for details.

    Input:
    R - A complex vector representing the X + iY positions of each particle on the surface 
    period - The periodicity of the wave in physical units

    Output:
    Ω - conformally mapped closed contour of the fluid surface
    =#

    Ω = exp.(- im * R);

    return Ω
end

function conformalDepth(h::Real)
    #=
    conformalDepth is a function that takes the real depth h and transforms it to the more useful conformal value H, which tends to 0 at inifinite depth. 
    To avoid working with infinite values, the infinite depth assumption is taken when the user inputs 0 as the real h.
    =#
    if iszero(h)
        H = 0.0
    else
        H = exp(-2 * h)
    end

    return H 
end

#=
The following DDI1 and DDI2 functions are used for taking the first and second order tangential derivatives with respect to the particle label ξ.
They do this according to Dold's weighted coefficients method for Lagrangian polynomial interpolation.
They are implemented in the main function with respect to conformal methods such that no offset or halo boundary is needed to deal with periodic boundary conditions.
=#
function DDI1(Ω::AbstractVector{<:Number}, offset)
    #=
    DDI1 is a function that uses an 11-point finite difference stencil to estimate the first derivative of Ω with respect to particle label ξ.
    =#
    N = length(Ω)
    Ω_ξ = similar(Ω)

    for i in 1:N
        Ω_ξ[i] = COEFFICIENTS1[1]*(Ω[mod1(i+1, N)]+offset*fld(i,N) - (Ω[mod1(i-1, N)]+offset*fld(i-2,N))) +
        COEFFICIENTS1[2]*(Ω[mod1(i+2, N)]+offset*fld(i+1,N) - (Ω[mod1(i-2, N)]+offset*fld(i-3,N))) +
        COEFFICIENTS1[3]*(Ω[mod1(i+3, N)]+offset*fld(i+2,N) - (Ω[mod1(i-3, N)]+offset*fld(i-4,N))) +
        COEFFICIENTS1[4]*(Ω[mod1(i+4, N)]+offset*fld(i+3,N) - (Ω[mod1(i-4, N)]+offset*fld(i-5,N))) +
        COEFFICIENTS1[5]*(Ω[mod1(i+5, N)]+offset*fld(i+4,N) - (Ω[mod1(i-5, N)]+offset*fld(i-6,N)))
    end   
    return Ω_ξ
    # First, extend Ω with 5 points on either side for periodicity considerations
    # Ωext = Vector{eltype(Ω)}(undef,N+10)
    # @inbounds Ωext[1:5] = Ω[N-4:N] .- offset 
    # @inbounds Ωext[6:N+5] = Ω
    # @inbounds Ωext[N+6:end] = Ω[1:5] .+ offset

    # Ω_ξ = similar(Ω)

    # @inbounds for i in 6:N+5
    #     Ω_ξ[i-5] = COEFFICIENTS1[1]*(Ωext[i+1] - Ωext[i-1]) + 
    #              COEFFICIENTS1[2]*(Ωext[i+2] - Ωext[i-2]) +
    #              COEFFICIENTS1[3]*(Ωext[i+3] - Ωext[i-3]) +
    #              COEFFICIENTS1[4]*(Ωext[i+4] - Ωext[i-4]) +
    #              COEFFICIENTS1[5]*(Ωext[i+5] - Ωext[i-5])
    # end

    # return Ω_ξ
end

function DDI2(Ω::AbstractVector{<:Number},offset)
    #=
    DDI2 is a function that uses an 11-point finite difference stencil to estimate the first derivative of Ω with respect to particle label ξ
    =#
    N = length(Ω)
    Ω_ξξ = similar(Ω)

    for i in 1:N
        Ω_ξξ[i] = COEFFICIENTS2[1]*Ω[i] +
        COEFFICIENTS2[2]*(Ω[mod1(i+1, N)]+offset*fld(i,N) + Ω[mod1(i-1, N)]+offset*fld(i-2,N)) +
        COEFFICIENTS2[3]*(Ω[mod1(i+2, N)]+offset*fld(i+1,N) + Ω[mod1(i-2, N)]+offset*fld(i-3,N)) +
        COEFFICIENTS2[4]*(Ω[mod1(i+3, N)]+offset*fld(i+2,N) + (Ω[mod1(i-3, N)]+offset*fld(i-4,N))) +
        COEFFICIENTS2[5]*(Ω[mod1(i+4, N)]+offset*fld(i+3,N) + (Ω[mod1(i-4, N)]+offset*fld(i-5,N))) +
        COEFFICIENTS2[6]*(Ω[mod1(i+5, N)]+offset*fld(i+4,N) + (Ω[mod1(i-5, N)]+offset*fld(i-6,N)))
    end   
    # Ω_ξξ = similar(Ω)
    # for i in 1:N
    #     points = [Ω[i],Ω[mod1(i+1, N)]+offset*fld(i,N) + Ω[mod1(i-1, N)]+offset*fld(i-2,N), Ω[mod1(i+2, N)]+offset*fld(i+1,N) + Ω[mod1(i-2, N)]+offset*fld(i-3,N), Ω[mod1(i+3, N)]+offset*fld(i+2,N) + Ω[mod1(i-3, N)]+offset*fld(i-4,N), Ω[mod1(i+4, N)]+offset*fld(i+3,N) + Ω[mod1(i-4, N)]+offset*fld(i-5,N), Ω[mod1(i+5, N)]+offset*fld(i+4,N) + Ω[mod1(i-5, N)]+offset*fld(i-6,N)]
    #     Ω_ξξ[i] = sum(COEFFICIENTS2 .* points)
    # end

    return Ω_ξξ
end

function ABMatrices(Ω::AbstractVector{<:Number}, Ω_ξ::AbstractVector{<:Number}, Ω_ξξ::AbstractVector{<:Number}, H::Real=0.0)
    #=
    The ABMatrices function sets up the A and B matrices from Dold eq. 4.13, using a necesary conditional double for loop.
    While it is a computationally expensive function, it is separated from the matrix inversion step because it is only needed once per timestep.

    Input:
    Ω - complex vector of conformally mapped positions of particles
    Ω_ξ - first-order derivative of Ω with respect to particle labels ξ
    N - number of particles
    H - parameter specified in Dold for bottom depth boundary condition. Defaulted to 0 for infinite depth.

    Output:
    A - Matrix
    B - Matrix
    ℵ - matrix of π I - B for use in NormalInversion
    =#
    N = length(Ω)
    if iszero(H)
        # create matrix of differences
        ΔΩ = Ω .- permutedims(Ω)
        # Compute off-diagonal elements first (diagonals will be Inf)
        C = Ω_ξ ./ (ΔΩ)
        # Fill in diagonal elements 
        C[diagind(C)] .= Ω_ξξ ./ (2 * Ω_ξ)
    else #
        C = zeros(Complex, N, N)    # C = A + iB
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

    ℵ = factorize(π*I - B)  # Identity matrix I from LinearAlgebra to subtract B from pi diagonal.

    return A,B,ℵ
end

function NormalInversion(ϕ::AbstractVector{<:Real}, A::AbstractMatrix{<:Real}, ℵ)
    #= 
    NormalInversion is a function that implements the matrix inversion method from Dold eq 4.13 in order to compute the normal derivative of the scalar velocity potential. The method is based on using a conformal mapping and the Cauchy integral theorem for solving the Laplacian equation. The subtlety lies in the issue that, for solely surface particles at b=0, there is no Lagrangian normal derivative as there are no particles above or below. For efficiency, due to the need of the tangential derivative, it is first computed and returned along with the normal derivative here.
        
    The implementation here is mathematically simplified from Dold's formula to take the form of an Ax = b matrix equation, and is thus quite compact and optimized.

    Input:
    ϕ - real vector of scalar velocity potential
    A - matrix built from ABMatrices
    ℵ - 

    Output:
    ϕ_ξ - real vector of tangential partial derivative of with respect to particle label ξ
    ϕ_ν - real vector of normal partial derivative scaled by
    =#

    ϕ_ξ = DDI1(ϕ,0)
    ϕ_ξξ = DDI2(ϕ,0)

    # Important: here the * is not element wise to get the sum A*ϕ_ξ for each one-element row entry of the resulting column vector, while the difference is element wise to subtract ϕ_ξξ[i] from each of the summed entries.
    b = ((A * ϕ_ξ) .- ϕ_ξξ)

    # Ax = b using the efficient \ operator, where x is the vector of tangential derivatives
    ϕ_ν = ℵ \ b

    return ϕ_ξ, ϕ_ν
end

function PhiTimeDer(R_ξ, ϕ_ξ, ϕ_ν, Y, p::SimulationParameters)
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
    ϕ_D = 0.5 .* (ϕ_ξ.^2 .+ ϕ_ν.^2) ./ abs.(R_ξ).^2 .- p.g .* Y
    ϕ_t = -0.5 .* (ϕ_ξ.^2 .+ ϕ_ν.^2) ./ abs.(R_ξ).^2 .- p.g .* Y

    return ϕ_D, ϕ_t
end

function TimeDerivatives(R_ξ::AbstractVector{<:Complex}, ϕ_x::AbstractVector{<:Real}, ϕ_y::AbstractVector{<:Real}, A::AbstractMatrix{<:Real}, ℵ,Y::AbstractVector{<:Real},p::SimulationParameters)
    #=
    TimeDerivatives is a function that computes up to the third Lagrangian time derivative
    of both the positions (x,y) and velocity potential ϕ
    =#

    # First derivatives of x,y just found from ϕ_x, ϕ_y which we already have
    # Evolution of ϕ comes from Bernoulli's equation at the free surface    
    DϕDt = 0.5 .* (ϕ_x.^2 .+ ϕ_y.^2) .- p.g.*Y

    # Get Eulerian time derivatives of velocities
    ϕξt, ϕνt = NormalInversion(-p.g.*Y .- 0.5 .*(ϕ_x.^2 .+ ϕ_y.^2), A, ℵ)
    ut, vt = RealPhi(R_ξ,ϕξt,ϕνt)

    # Next, get gradient of velocity fields
    xξ = real(R_ξ)
    yξ = imag(R_ξ)
    spacing = xξ.^2 .+ yξ.^2
    uξ = DDI1(ϕ_x,0)
    vξ = DDI1(ϕ_y,0)
    utξ = DDI1(ut,0)
    vtξ = DDI1(vt,0)
    ux = (uξ.*xξ .- vξ.*yξ)./spacing
    vy = -ux
    vx = (uξ.*yξ .+ vξ.*xξ)./spacing
    uy = vx
    utx = (utξ.*xξ .- vtξ.*yξ)./spacing
    vty = -utx
    vtx = (utξ.*yξ .+ vtξ.*xξ)./spacing
    uty = vtx
    uxξ = DDI1(ux,0)
    vxξ = DDI1(vx,0)
    uxx = (uxξ.*xξ .- vxξ.*yξ)./spacing
    vxx = (uxξ.*yξ .+ vxξ.*xξ)./spacing
    

    # Next, use this to get Lagrangian time derivatives
    DuDt = ut .+ ϕ_x .* ux .+ ϕ_y .* uy
    DvDt = vt .+ ϕ_x .* vx .+ ϕ_y .* vy
    D2ϕDt2 = ϕ_x .* DuDt .+ ϕ_y .* DvDt .- p.g .* ϕ_y

    # Now do third order, need to first get utt Eulerian
    ϕttξ, ϕttν = NormalInversion(-ϕ_x .* ut .- ϕ_y .* vt,A, ℵ)
    utt,vtt = RealPhi(R_ξ,ϕttξ,ϕttν)

    # Compute third order time derivatives
    D2uDt2 = utt .+ 2 .*(ϕ_x .* utx .+ ϕ_y.*uty) .+ ut.*ux .+ vt.*vx .+ (ux.^2 .+ vx.^2).*ϕ_x .+ (ϕ_x.^2 .- ϕ_y.^2).* uxx .+ 2 .* ϕ_x .* ϕ_y .* vxx
    D2vDt2 = vtt .+ 2 .*(ϕ_x .* vtx .+ ϕ_y.*vty) .+ ut.*vx .- vt.*ux .+ (ux.^2 .+ vx.^2).*ϕ_y .+ (ϕ_x.^2 .- ϕ_y.^2).* vxx .- 2 .* ϕ_x .* ϕ_y .* uxx
    D3ϕDt3 = DuDt.^2 .+ DvDt.^2 .+ ϕ_x .* D2uDt2 .+ ϕ_y .* D2vDt2 .- p.g.*DvDt

    return DϕDt, DuDt, DvDt, D2ϕDt2, D2uDt2, D2vDt2, D3ϕDt3
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

function LagrangeInterpolant(ϕ::AbstractVector{<:Real},tReal::AbstractVector{<:Real},l::Int)
    #=
    LagrangeInterpolant is a function which, given l previous points in time,
    forms a Lagrange interpolating polynomial of order l. In practice, this is
    performed on the third order time derivatives of X,Y,ϕ. This way, we can 
    use the coefficients to predict up to the (l+3)th order time derivatives
    using this interpolating polynomial in time. We then output those higher
    order derivatives to use for Taylor timestepping.

    Inputs:
    ϕ - a vector of of function values to interpolate 
    tReal - a vector of times where each function value is taken
    =#
    # Shift t-values so that t[end] = 0 (makes derivatives simpler to compute)
    t = tReal .- tReal[end]
    # Create empty vector of derivative estimates (up to 5th)
    dϕ = zeros(5)
    n = l + 1

    # First, compute weights of Lagrange interpolating polynomial 
    w = ones(n)
    for i ∈ 1:n 
        for j ∈ 1:n 
            i == j && continue
            w[i] /= (t[i]-t[j])
        end 
    end 

    # Compute higher order derivatives given weights 
    for i ∈ 1:n 
        δ = -t[i]
        dϕ[1] += ϕ[i] * w[i]
        pow = 1.0 
        for m ∈ 2:n 
            pow *= δ
            dϕ[m] += ϕ[i] * w[i] * pow 
        end 
    end 

    for m ∈ 1:n 
        dϕ[m] *= factorial(m)
    end 

    return dϕ


    # if l == 0 # i.e. first timestep 
    #     nothing # if only one point, can't estimate further
    # end

    # if l == 1  #i.e. two point linear fit 
    #     p1 = LagrangeInterpolation(t,ϕ)
    #     dϕ[1] = derivative(p1,1)(0)
    # end 

    # if l == 2  #i.e. three point quadratic fit 
    #     p2 = LagrangeInterpolation(t,ϕ)
    #     dϕ[1] = derivative(p2,1)(0)
    #     dϕ[2] = derivative(p2,2)(0)
    # end

    # if l == 3 #i.e. four point cubic fit 
    #     p3 = LagrangeInterpolation(t,ϕ)
    #     dϕ[1] = derivative(p3,1)(0)
    #     dϕ[2] = derivative(p3,2)(0)
    #     dϕ[3] = derivative(p3,3)(0)
    # end

    # if l == 4 #i.e. five point quartic fit 
    #     p4 = LagrangeInterpolation(t,ϕ)
    #     dϕ[1] = derivative(p4,1)(0)
    #     dϕ[2] = derivative(p4,2)(0)
    #     dϕ[3] = derivative(p4,3)(0)
    #     dϕ[4] = derivative(p4,4)(0)
    # end

    # if l == 5 #i.e. six point quintic fit 
    #     p5 = LagrangeInterpolation(t,ϕ)
    #     dϕ[1] = derivative(p5,1)(0)
    #     dϕ[2] = derivative(p5,2)(0)
    #     dϕ[3] = derivative(p5,3)(0)
    #     dϕ[4] = derivative(p5,4)(0)
    #     dϕ[5] = derivative(p5,5)(0)
    # end

    # return dϕ
end

function LagrangeInterpolation(t,ϕ)
    n = length(t)
    Lagrange = [Polynomial(1.0) for i=1:n]
    for i in 1:n
        for j in 1:n
            if i != j
                Lagrange[i] *= Polynomial([-t[j], 1]) / (t[i] - t[j])
            end
        end
    end
    return sum(ϕ[i] * Lagrange[i] for i in 1:n)
end



function RealPhi(R_ξ::AbstractVector{<:Complex}, ϕ_ξ::AbstractVector{<:Real}, ϕ_ν::AbstractVector{<:Real})
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

function smooth(Ω::AbstractVector{<:Number},offset)
    N = length(Ω)
    Ω_sm = similar(Ω)

    for i in 1:N
        
        points = [Ω[i],Ω[mod1(i+1, N)]+offset*fld(i,N) + Ω[mod1(i-1, N)]+offset*fld(i-2,N), Ω[mod1(i+2, N)]+offset*fld(i+1,N) + Ω[mod1(i-2, N)]+offset*fld(i-3,N), Ω[mod1(i+3, N)]+offset*fld(i+2,N) + Ω[mod1(i-3, N)]+offset*fld(i-4,N), Ω[mod1(i+4, N)]+offset*fld(i+3,N) + Ω[mod1(i-4, N)]+offset*fld(i-5,N), Ω[mod1(i+5, N)]+offset*fld(i+4,N) + Ω[mod1(i-5, N)]+offset*fld(i-6,N), Ω[mod1(i+6, N)]+offset*fld(i+5,N) + Ω[mod1(i-6, N)]+offset*fld(i-7,N), Ω[mod1(i+7, N)]+offset*fld(i+6,N) + Ω[mod1(i-7, N)]+offset*fld(i-8,N)]

        δM = SMOOTHCOEFFICIENTS' * points
        Ω_sm[i] = Ω[i] - δM
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


function computeEnergy(sol,n,L=2π)
    #treat initial case differently 
    # x = sol(0)[1:n]
    # y = sol(0)[n+1:2*n]
    # ϕ = sol(0)[2*n+1:end]
    # xξ = DDI1(x,n,L,1)
    # yξ = DDI1(y,n,0,1)
    # du = zeros(length(sol(0)))
    # TimeStep(du,sol(0),[n,L,0,false],0)
    # ẋ = du[1:n]
    # ẏ = du[n+1:2*n]
    # B_ξ = ẋ.*yξ .- ẏ.*xξ
    # integrand = -1/2 * ϕ .* B_ξ
    # energy = [sum(integrand)*(2π)/n + sum(GRAVITY/2 * (y).^2 .* xξ)*2π/n]
    # MWL_check = [sum(y.*xξ)/L]
    # phasespd = [median(B_ξ./yξ)];
    KE = [];
    PE = [];
    MWL_check = [];
    phasespd = [];
    for t ∈ sol.t[1:end]
        x = sol(t)[1:n]
        y = sol(t)[n+1:2*n]
        ẋ = sol(t,Val{1})[1:n]
        ẏ = sol(t,Val{1})[n+1:2*n]
        ϕ = sol(t)[2*n+1:end]
        xξ = DDI1(x,L)
        yξ = DDI1(y,0)
        #R = x .+ im.*y
        #Ω = conformalMap(R)
        #Ω_ξ = DDI1(Ω, n,0)
        #R_ξ = (im ./ Ω) .* Ω_ξ
        #q = ẋ .+ im.*ẏ
        #ϕ_n = imag.(q .* abs.(R_ξ) ./ R_ξ)
        # Other way of getting KE from Balk
        MWL = sum(y.*xξ)/ n
        B_ξ = ẋ.*yξ .- ẏ.*xξ
        integrand = -1/2 * ϕ .* B_ξ
        #KE = simpsons_rule_periodic(x,ϕ.* ϕ_n/2)
        Kinetic = sum(integrand)*(2π)/n
        #println(KE)
        Potential = sum(GRAVITY/2 * (y.-MWL).^2 .* xξ)*2π/n
        append!(MWL_check,MWL) # Eulerian MWL
        #append!(MWL_check,sum(y)/n) # Lagrangian MWL
        append!(KE,Kinetic)
        append!(PE,Potential)
        append!(phasespd,median(B_ξ./yξ))
        #append!(momentum, p)
    end
    return KE, PE , MWL_check, phasespd
end

function computeEnergyDold(xvals,yvals,ϕvals,time,N,c,L=2π)
    KE = [];
    PE = [];
    MWL_check = [];
    phasespd = [];
    momentum = [];
    for t ∈ enumerate(time)
        x = xvals[t[1],:]
        y = yvals[t[1],:]
        ϕ = ϕvals[t[1],:]
        xξ = DDI1(x,L)
        yξ = DDI1(y,0)
        (ẋ, ẏ, _, _, _, _, _, _, _) = fixedTimeOperations(N, x, y, ϕ, L, 0.0)
        
        MWL = sum(y.*xξ)/ N
        B_ξ = ẋ.*yξ .- ẏ.*xξ
        integrand = -1/2 * ϕ .* B_ξ
        Kinetic = sum(integrand)*(2π)/N
        Potential = sum(GRAVITY/2 * (y).^2 .* xξ)*2π/N
        momIntegral = sum(ẋ.*xξ)/N 
        append!(MWL_check,MWL)
        append!(KE,Kinetic)
        append!(PE,Potential)
        append!(phasespd,median(B_ξ./yξ))
        append!(momentum,momIntegral)
    end
    return KE, PE, MWL_check, phasespd, momentum
end
