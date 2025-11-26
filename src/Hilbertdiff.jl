using DrWatson
@quickactivate "Castawave"
using LaTeXStrings
using Printf
using MAT
using Trapz
using Plots
using JLD2
plotlyjs()

include(projectdir()*"/src/MainSolver.jl")
include(projectdir()*"/src/ClamondIC.jl")
include(projectdir()*"/src/DoldMono.jl")

N = 128;
A = 0.3
X,Y,ϕ,c = getIC(Inf,A,N÷2)
tf = 4.0;
Δt = 1e-4
L = 2π
h = 0.0
tol = 1e-6
smoothed = true

Xfull, Yfull, ϕfull, t = @time runSim(N, X, Y, ϕ, Δt, Float64(tf),L,h,ϵ = tol,smoothing=smoothed)

anim = @animate for i ∈ enumerate(t)
    ϕ = ϕfull[i[1],:]
    X = Xfull[i[1],:]
    Y = Yfull[i[1],:]
    xξ = DDI1(X,N,L,1)
    yξ = DDI1(Y,N,0,1)
    (ẋ, ẏ, _, ẍ, ÿ, _, _, _, _) = fixedTimeOperations(N, X, Y, ϕ, L, 0.0);
    ψξ = ẋ.*yξ .- ẏ.*xξ
    ϕξ = DDI1(ϕ,N,0,1)

    #r = abs.(ẍ .+ im * ÿ .+ im*GRAVITY) ./ abs.(xξ .+ im*yξ)
    # Ω = conformalMap(X .+ im*Y)
    # R_ξ = DDI1(X,N,L,1) .+ im*DDI1(Y,N,0,1)
    # R_ξξ = DDI2(X,N,L,1) .+ im*DDI2(Y,N,0,1)
    # Ω_ξ = DDI1(Ω,N,0,0)
    # Ω_ξξ = DDI2(Ω,N,0,0)

    # H = conformalDepth(h)

    # # The matrix method described in Dold is used to find the normal derivative of the potential.
    # A, B, ℵ = ABMatrices(Ω, Ω_ξ, Ω_ξξ, N, H)
    # ϕ_ξ, ϕ_ν = NormalInversion(ϕhalf, A, ℵ, N)
    # ψξξ = DDI1(ψξ,N,0,1);
    # b = A * ψξ .- ψξξ;
    # ψν = ℵ \ b
    diff = imag.(hilbert(ψξ)) .- ϕξ
    diff2 = imag.(hilbert(ϕξ)) .+ ψξ
    plot(X,diff,label = "H(ψξ) - ϕξ",title = @sprintf("Hilbert deviations Time: %.1f s", i[2]),ylims=(-0.03,0.03),xlabel="x (m)")
    plot!(X,diff2,label = "H(ϕξ) + ψξ")
    #plot(r,label="r")

end every 10
gif(anim, projectdir()*"/plots/RK4noSmooth.gif", fps=10)

plot(Xfull[100,:],Yfull[100,:])
z = Xfull[1,:] .+ im*Yfull[1,:]
zξ = DDI1(Xfull[1,:],N,2π,1) .+ im*DDI1(Yfull[1,:],N,0,1)
plot(real(z),imag(z))
Ω = conformalMap(z)
Ωξ = conformalMap(zξ)
plot(real(Ω),imag(Ω))
scatter(real(Ωξ),imag(Ωξ))
using Trapz
ξ = 1:N
ξ_p = collect(1:N) .-0.2*im 
zint = ComplexF64.(zeros(N))
for i ∈ 1:N
    zint[i] = -1/(2π * im) * (trapz(ξ ,Ωξ ./(ξ .- (ξ_p[i]) )))
end
scatter(real(zξ),imag(zξ))
scatter!(real(zint),imag(zint))
scatter(real(Ω),imag(Ω))
scatter!(real(zint),imag(zint))
z = im*log.(Ωξ)
zintObs = im*log.(zint)
scatter(real(zξ),imag(zξ))
scatter!(real(zintObs),imag(zintObs))
# Try doing this in circular domain 
half = length(t)÷ 2
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

plot(Ω,aspect_ratio=1)
plot(real.(Ω_ξ),imag.(Ω_ξ),aspect_ratio=1)
plot(abs.(fft(real.(Ω))).^2)
plot(abs.(fft(imag.(Ω))).^2)

plot(imag.(Ω).+imag.(hilbert(real.(Ω))))
plot!(-imag.(hilbert(real.(Ω))))

plot(imag.(hilbert(ϕξ)) .+ ψξ,label="H(ϕξ) + ψξ")
plot(imag.(hilbert(xξξ)).+yξξ,label= "H(xξξ) + yξξ")
plot(imag.(hilbert(yξξ)).-xξξ)

plot!(imag.(hilbert(ψξ)) .- ϕξ,label="H(ψξ) - ϕξ")
plot!(imag.(hilbert(yξ)).-xξ,label="H(yξ).- xξ")
ds = xξ.^2 .+ yξ.^2 
dsξ = DDI1(ds,N,0,1)
plot(dsξ)
plot(ds)

plot(ψξ)
plot!(c*yξ)

Y_ξ, Y_ν = NormalInversion(Y, A, ℵ, N)
plot(Y_ν.-ϕξ/c)
plot(imag.(hilbert(Y_ν)).-xξ)
plot!(xξ)
plot(c*Y_ν .- ψν)

X_ξ = DDI1(X, N,L,1)
X_ξξ = DDI2(X, N,L,1)

# Important: here the * is not element wise to get the sum A*ϕ_ξ for each one-element row entry of the resulting column vector, while the difference is element wise to subtract ϕ_ξξ[i] from each of the summed entries.
b = ((A * X_ξ) .- X_ξξ)

# Ax = b using the efficient \ operator, where x is the vector of tangential derivatives
X_ν = ℵ \ b
plot(X_ν,label="Xν")
plot!(X_ξ,label="Xξ")
plot!(Y_ξ,label="Yξ")
plot!(Y_ν,label="Yν")

plot(imag.(hilbert(ϕξ)))
plot!(-ψξ)
plot(imag.(hilbert(ϕξ)) + ψξ)
plot!(imag.(hilbert(xξ)).+yξ)

plot(ϕν)
plot!(ϕ)

# Look at curvature
yξξ = DDI2(Y,N,0,1)
xξξ = DDI2(X,N,2π,1)
κ = (yξξ.*xξ .- xξξ.*yξ)./(ds.^(3/2))
plot(κ)
plot!(sqrt.(ds).*10)
v2 = ẋ.^2 + ẏ.^2
plot(v2)
plot!(κ)
plot((ẍ.*(-yξ).+ ÿ.*xξ)./sqrt.(ds) .- κ)
plot((ẍ.*(-yξ).+ (ÿ).*xξ)./sqrt.(ds),label="Normal acceleration")
plot(-Y)
plot!(sqrt.(ds).*10)
plot(ẍ.*xξ .+ (ÿ .+ GRAVITY).*yξ)
plot!(-1/2*c^2*Y)
plot!(κ.*ds)

# Try new way to get ψ
G(z,z0) = 1/(2π)* log.(abs.(z .- z0))
∇G(z,z0) = 1/(2π)./ conj.(z .- z0)
n̂(xξ,yξ) = (-yξ .+ im.*xξ)./(sqrt.(xξ.^2 .+ yξ.^2))
ŝ(xξ,yξ) = (xξ .+ im*yξ)./(sqrt.(xξ.^2 .+ yξ.^2))

Z = X .+ im.*Y

real.(∇G(Z[1],Z[2]).* conj.(ŝ(xξ[1],yξ[2]))).* ϕ[1]


Ω = conformalMap(X .+ im*Y)
R_ξ = DDI1(X,N,L,1) .+ im*DDI1(Y,N,0,1)
R_ξξ = DDI2(X,N,L,1) .+ im*DDI2(Y,N,0,1)
Ω_ξ = DDI1(Ω,N,0,0)
Ω_ξξ = DDI2(Ω,N,0,0)

H = conformalDepth(h)

# The matrix method described in Dold is used to find the normal derivative of the potential.
A, B, ℵ = ABMatrices(Ω, Ω_ξ, Ω_ξξ, N, H)
ϕ_ξ, ϕ_ν = NormalInversion(ϕ, A, ℵ, N)
ψξξ = DDI1(ψξ,N,0,1);
b = A * ψξ .- ψξξ;
ψν = ℵ \ b

b = -A * ϕ .+ DDI1(ϕ,N,0,1)
ψ = (2* I*π - B) \ b

plot(ψξ)
plot!(DDI1(ψ,N,0,1))
plot(ψξ .+ DDI1(ψ,N,0,1))
plot!(imag.(hilbert(ϕξ)).+DDI1(ψ,N,0,1))
plot(-ψ)

b = A * ϕξ .- DDI2(ϕ,N,0,1)
ϕν = ℵ \ b

plot(ϕν)
plot!(-ψξ)
plot!(ϕξ)
plot!(ψν)