using DrWatson
@quickactivate "Castawave"

using LinearAlgebra 
using JLD2
using DifferentialEquations
using FFTW
using ForwardDiff

include(projectdir()*"/src/Constants.jl")
include(projectdir()*"/src/HelperFunctions.jl")

function FFT_surface(Y)
    #=
    This function takes the Fourier transform of Y(ξ) at a fixed time and
    outputs Yk values, the Fourier coefficients at each mode.

    =#
    
    # first, take the fourier transform, shifting so that it goes from (-K,K)
    Yfft = fftshift(fft(Y))
    K = length(Y)÷2
    freq = range(-K,K)
    # Deal with nyquist frequency, split evenly into -K and K modes 
    Yfft[1] /= 2
    push!(Yfft,Yfft[1])
    
    return Yfft

end

f(m,k) = sign(k) == sign(m-k) ? 0 : 1

bracket(j) = j == 0 ? 1 : abs(j)

P(m,n,Yk) = 2/sqrt(bracket(m)) * (f(m,n)*abs(m-n)*Yk[m-n] .- abs(m*n)*Yk[m]*Yk[-n])


####### Start with one mode 
Att(A,B,At,Bt) = 2*(-8 * (-1 + At^2)*B + A*(1 + 8*(-8 + 7*At^2)*B^2 + 16*At*Bt) - 4*A^2*(B + 48*(-1 + At^2)*B^3 + 12*At*B*Bt) + 8*A^3*(32*(-1+At^2)*B^4 - Bt^2))/((1 - 4*A*B)^2)
Btt(A,B,At,Bt) = -2*B*(-1+8*(-4+5*At^2)*B^2 - 16*At*Bt + 4*A*B*(1-48*(-1+At^2)*B^2 + 12*At*Bt) + 8*A^2*(32*(-1+At^2)*B^4 + Bt^2))/((1-4*A*B)^2)

function odes!(du,u,p,t)
    A,At,B,Bt = u
    du[1] = At
    du[2] = Att(A,B,At,Bt)
    du[3] = Bt
    du[4] = Btt(A,B,At,Bt)
end

# generating initial conditions 
a0 = 0.001;
c = sqrt(1-4*a0^2);
v0 = im*c*a0

u0 = [a0,v0,conj(a0),conj(v0)]
tspan=(0.0,4.0)


prob = ODEProblem(odes!,u0,tspan)
sol =  solve(prob,Vern9(), abstol = 1e-16, reltol = 1e-16)

sols = hcat(sol.u...)
plot(sol.t,abs.(sols[1,:]),ylabel="|Y₁|",xlabel="t",legend=false)

A,At,B,Bt = u0
Att(A,B,At,Bt)
Btt(A,B,At,Bt)

Att(0.1,0.0,0.2,0)
Btt(0.1,0.0,0.2,0.0)