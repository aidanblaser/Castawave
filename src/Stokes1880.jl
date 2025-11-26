using DrWatson
@quickactivate "Castawave"
using Plots
using Polynomials

include(projectdir()*"/src/ClamondIC.jl")

# Get initial conditions from Clamond
N = 512;
A = 0.2;
X,Y,ϕ,c = getIC(Inf,A,N÷2)
ϕsub = ϕ[97:420]
Xsub = X[97:420]
Ysub = Y[97:420]
g = 9.81;
# Rescale ϕ to 2π domain 
ϕp = (ϕ*π/maximum(ϕ))
plot(ϕp)

# Define surface in terms of ϕ from Stokes (1880) eqns (17) + (18)
x(ϕ,b) = -ϕ .+ b * sin.(ϕ) .- (b^2 + 1/2*b^4)*sin.(2 .*ϕ) .+ (3/2*b^3 + 19/12*b^5)*sin.(3 .* ϕ) .- 
8/3 * b^4 *sin.(4 .* ϕ) .+ 125/24 * b^5 * sin.(5 .* ϕ)

y(ϕ,b) = b*cos.(ϕ) .- (b^2 + 1/2 * b^4)*cos.(2 .* ϕ) .+ (3/2*b^3 + 19/12 * b^5)*cos.(3 .* ϕ) .-
8/3 * b^4 * cos.(4 .* ϕ) .+ 125/24 * b^5 * cos.(5 .* ϕ)

# Determine b from max height (solve quintic polynomial, Galois is not pleased)
p = Polynomial([-A,1,0,3/2,0,163/24]);
r = roots(p);
# Only take real root 
b = Float64(r[3])

# Put into x and y (taking into account phase shifts)
y0 = 1/2*A^2;
xsub = x(ϕp,b);
ysub = y(ϕp,b).+y0;

plotlyjs()
scatter(xsub,ysub,xlabel="x (m)",ylabel="y (m)",title="Free Surface",label="5th order Stokes solution")
scatter!(X,Y,label="Clamond")


plot(xsub)
plot(ysub)


scatter(X,Y)
scatter!(X,ϕ)
plot(ϕ)

