using DrWatson
@quickactivate "Castawave"
using Test
using MAT
using LaTeXStrings
using Printf
using Sundials

#= Test 1: Does the phase speed match theory?

For this test, we will be simulating waves of varying steepness up
until Ak = 0.42. By using maximum tracking, we will determine the
phase speed, and plot it as a function of steepness.

=# 

include(projectdir()*"/src/MainSolver.jl")
include(projectdir()*"/src/ClamondIC.jl")

# Initial Parameters 
# Range of amplitudes
Arange = range(start=0.001,stop=0.42,length=10);
# Set wavenumber to 1, won't change scaling
k = 1;
L = 2π
# Number of points
N = 128;
# Set smoothing
smoothed = true
# Timestep and final time
Δt = 0.001;
tf = Float64(2000);
alg = ImplicitEulerExtrapolation(min_order=5,autodiff=false)
alg = AutoVern9(ESDIRK547L2SA2(autodiff=false))
alg = AutoVern9(Rodas4P(autodiff=false))
alg = CVODE_BDF()
alg = Rodas4P(autodiff=false)
alg = Kvaerno4(autodiff=false)
alg = ImplicitEuler(autodiff=false)
alg = Vern7()
alg = CVODE_Adams()

A = 0.2;
X = [(α * L / N) .- A*sin(L*α/N) .- A^3*k^2*sin(L*α/N) .- A^4*k^3 / 3 * sin(2*L*α/N) for α in 1:N];
Y = [(cos(L * α / N )) * A + 1/6*A^4*k^3*cos(L*α/N) .+ (A^2*k / 2).+ A^4*k^3 * 1/2 for α in 1:N];
ϕ = [sqrt(9.81/k) * A  * exp.(k*Y[α]) * sin(k*X[α]) for α in 1:N];
Y .-= MWL


using MAT
vars = matread(projectdir()*"/data/ClamondA2N128.mat")
X = vec(vars["X"])
Y = vec(vars["Y"]);
ϕ = vec(vars["F"])*sqrt(9.81);
MWL = simpsons_rule_periodic(X,Y)
Y .-=MWL

# Interpolate onto more even grid
using CubicSplines
splineY = CubicSpline(X,Y)
splineϕ = CubicSpline(X,ϕ)
X = collect(range(minimum(X),stop=maximum(X),length=N))
Y = splineY[X]
ϕ = splineϕ[X]

using DSP
sol = runSim(N,X,Y,ϕ,Δt,tf,L,smoothing=smoothed,alg = alg,relabelling=true)


# Get X, Y, ϕ values

t = 0:10*Δt:tf
xvals = zeros(length(t),N);
yvals = zeros(length(t),N);
ϕvals = zeros(length(t),N);
for i ∈ 1:length(t)
    xvals[i,:] = sol(t[i])[1:N]
    yvals[i,:] = sol(t[i])[N+1:2*N]
    ϕvals[i,:] = sol(t[i])[2*N+1:3*N]
end

scatter(sol(2000)[1:N],sol(2000)[N+1:2*N])




sol(10)[1:N]
# Interpolate
using CubicSplines
maxVec = []
for i ∈ t
    xval = sol(i)[1:N]
    yval = sol(i)[N+1:2*N]
    spline = CubicSpline(xval,yval)
    xrange = range(minimum(xval),stop=maximum(xval),length=10000)
    η = spline[xrange]
    push!(maxVec,xrange[argmax(η)])
end

c = diff(maxVec)./Δt 
cfilt = max.(c,2.7) 
sum(cfilt)./ length(cfilt)
sqrt(9.81)*(1 + 1/2 *(0.1^2) + 1/8 * (0.1)^4)
plot(t[1:end-1],cfilt,ylims=(0,5))
plot!([0,10],[sqrt(9.81)*(1 + 1/2 *(0.3^2) + 1/8 * 0.3^4),sqrt(9.81)*(1 + 1/2 *(0.3^2)+ 1/8 * 0.3^4)])


gr()
function visualize(interval::Int, fps::Int)
    anim = @animate for i ∈ t
        xvals = mod.(sol(i)[1:N],2π)
        yvals = sol(i)[N+1:2*N]
        ϕvals = sol(i)[2*N+1:3*N]

        scatter(xvals,yvals, legend = false,
         framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1,
          dpi = 300, xlabel=L"x \,(m)",ylabel=L"z \,(m)", title= @sprintf("Time: %.1f s", i),
          xlims=(0,L),ylims = (-1.9,1.9),aspect_ratio=1)
        #plot!([maxVec[Int(t÷Δt + 1)],maxVec[Int(t÷Δt + 1)]],[-2,2],linewidth=3)
    end every interval
    gif(anim, projectdir()*"/plots/RK4noSmooth.gif", fps=fps)
end
visualize(20, 30)

energy, MWL_check, phasespd = computeEnergy(sol,N,Δt*100,tf)
gr()
plot(0:Δt*100:tf,energy .- energy[1],xlabel=L"t \, (s)",ylabel =L"E - E_0 \quad (J/kg)",legend=false,title= L"E_0: \quad "*@sprintf("%.4f ", energy[1])*L"(J/kg)")
plot(0:Δt*100:tf,MWL_check,xlabel=L"t \, (s)",ylabel ="MWL"*L" \quad (m)",legend=false)
plot(0:Δt:tf, momentum)

# Run simulations for range of steepnesses
for A ∈ Arange
    X = [(α * L / n) - 0*A*sin(k*α/N) - 0*A^3*k^2*sin(k*α/N) - 0*A^4*k^3 / 3 * sin(2*k*α/N) for α in 1:N]
    Y = [(cos(k * α / N )) * A + 0*1/6*A^4*k^3*cos(k*α/N) .+ 0*(A^2*k / 2).+ 0*A^4*k^3 * 1/2 for α in 1:N]
    ϕ = [sqrt(9.81/k) * A   * sin(k*X[α]) for α in 1:N]
    sol = runSim(N,X,Y,ϕ,Δt,tf,smoothing=smoothed)
end

ϕvals


t = 0
    x = sol(t)[1:n]
    y = sol(t)[n+1:2*n]
    ẋ = sol(t,Val{1})[1:n]
    ẏ = sol(t,Val{1})[n+1:2*n]
    ϕ = sol(t)[2*n+1:end]
    R = x .+ im.*y
    Ω = conformalMap(R)
    Ω_ξ = DDI1(Ω, n)
    R_ξ = (im ./ Ω) .* Ω_ξ
    q = ẋ .+ im.*ẏ
    ϕ_n = imag.(q .* abs.(R_ξ) ./ R_ξ)
    # Other way of getting KE from Balk 
    B_ξ = ẋ.*imag.(R_ξ) .- ẏ.*real.(R_ξ)
    ϕ_ξ = ẋ.*real.(R_ξ) .+ ẏ.*imag.(R_ξ)

using DSP
plot((hilbert(B_ξ)))
plot!(ϕ_ξ)


# Check eta stuff
using CubicSplines
ηvec = []
for i∈t
    x = mod.(sol(i)[1:N],2π)
    y = sol(i)[N+1:2*N]
    indices = sortperm(x)
    spline = CubicSpline(x[indices],y[indices])
    push!(ηvec,[spline,i])
end

ηvec[1]
ranges = range(0.08,2π-0.08,length=10000)
plot(ranges, ηvec[1000][1][ranges])

omega(A,k,g) = sqrt(g*k)*(1 + 1/2*A^2*k^2 + 1/8*A^4*k^4)
η(x,t,A,k,g) = A * cos.(k*x .- omega(A,k,g)*t)*(1 - 1/16*A^2*k^2) + (1/2*A^2*k + 5/6*A^4*k^3)*cos.(2*k*x .- 2*omega(A,k,g)*t) + 3/8*A^3*k^2*cos.(3*k*x .- 3*omega(A,k,g)*t) + 1/3*A^4*k^3*cos.(4*k*x .- 4*omega(A,k,g)*t)

# Comparison
for i ∈ 1:length(t) 
    derived = ηvec[i][1][ranges]
    theory = η(ranges,ηvec[i][2],0.2,1,9.81)
end

i = 10000
plotlyjs()
derived = ηvec[i][1][ranges];
theory = η(ranges,ηvec[i][2],0.1979,1,9.81);
plot(ranges,derived,label="Numerics")
plot!(ranges,theory,label="4th order theory")


ηvec[10000][2]

X,Y,ϕ,c =   getIC(Inf,0.2,64)
scatter(X,Y)
