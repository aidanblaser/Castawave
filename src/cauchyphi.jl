using DrWatson
@quickactivate "Castawave"
using Plots

include(projectdir()*"/src/MainSolver.jl")
include(projectdir()*"/src/ClamondIC.jl")
include(projectdir()*"/src/DoldMono.jl")

N = 1024;
A = 0.3
X,Y,ϕ,c = getIC(Inf,A,N÷2)
tf = 4.0;
Δt = 1e-4
L = 2π
h = 0.0
tol = 1e-6
smoothed = true

Xfull, Yfull, ϕfull, t = @time runSim(N, X, Y, ϕ, Δt, Float64(tf),L,h,ϵ = tol,smoothing=smoothed)

# Try doing this in circular domain 
half = length(t)÷ 2
half = 150
X = Xfull[half,:]
Y = Yfull[half,:]
ϕ = ϕfull[half,:]
Ω = conformalMap(X .+ im*Y)
R_ξ = DDI1(X,N,L,1) .+ im*DDI1(Y,N,0,1)
R_ξξ = DDI2(X,N,L,1) .+ im*DDI2(Y,N,0,1)
Ω_ξ = DDI1(Ω,N,0,0)
Ω_ξξ = DDI2(Ω,N,0,0)


xξ = DDI1(X,N,L,1)
yξ = DDI1(Y,N,0,1)
xξξ = DDI2(X,N,L,1)
yξξ = DDI2(Y,N,0,1)
(ẋ, ẏ, _, ẍ, ÿ, _, _, _, _) = fixedTimeOperations(N, X, Y, ϕ, L, 0.0);
ψξ = ẋ.*yξ .- ẏ.*xξ
ϕξ = DDI1(ϕ,N,0,1)
ds = sqrt.(xξ.^2 .+ yξ.^2 )

plot(ϕξ./ds)
plot!(ψξ./ds)

A, B ,ℵ = ABMatrices(Ω,Ω_ξ,Ω_ξξ,N,0.0)
ϕξ, ϕν = NormalInversion(ϕ,A,ℵ,N)
plot(ψξ)
plot!(-ϕν)
plot(ψξ .+ ϕν)

plot(abs.(Ω_ξ))

function ABMatricesAdj(Ω, Ω_ξ, Ω_ξξ, N, H=0.0)
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

    if iszero(H)
        for ξ_p in 1:N
            for ξ in 1:N
                if ξ_p == ξ     
                    C[ξ,ξ_p] = (Ω_ξξ[ξ]) / (2 * Ω_ξ[ξ])
                else
                    C[ξ,ξ_p] = (Ω_ξ[ξ_p]) / (Ω[ξ] - Ω[ξ_p])
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

    ℵ = factorize(π*I - B)  # Identity matrix I from LinearAlgebra to subtract B from pi diagonal.

    return A,B,ℵ
end

A, B, ℵ = ABMatricesAdj(Ω,Ω_ξ,Ω_ξξ,N,0.0)



function NormalInversionAdj(ϕ, A, ℵ, N)
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

    ϕ_ξ = DDI1(ϕ, N,0,1)
    ϕ_ξξ = DDI2(ϕ,N,0,1)

    # Important: here the * is not element wise to get the sum A*ϕ_ξ for each one-element row entry of the resulting column vector, while the difference is element wise to subtract ϕ_ξξ[i] from each of the summed entries.
    b = -(A * ϕ .- ϕ_ξ)

    # Ax = b using the efficient \ operator, where x is the vector of tangential derivatives
    ψ = ℵ \ b

    return ϕ, ψ
end

ϕ,ψ = NormalInversionAdj(ϕ,A,ℵ,N)

plot(ϕ)
plot!(ψ)


ψξest = DDI1(ψ,N,0,1)

plot(ψξ)
plot!(ψξest)
plot!(-imag.(hilbert(ϕξ)))

plot(ψξ .- ψξest)
plot(ds)