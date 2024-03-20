using DrWatson
@quickactivate "Castawave"
using Test
using MAT
using LaTeXStrings
using Printf

#= Test 1: Does the phase speed match theory?

For this test, we will be simulating waves of varying steepness up
until Ak = 0.42. By using maximum tracking, we will determine the
phase speed, and plot it as a function of steepness.

=# 

include(projectdir()*"/src/MainSolver.jl")

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
tf = Float64(1);
alg = ImplicitEulerExtrapolation(min_order=5,autodiff=false)
alg = Vern9()

A = 0.2;
X = [(α * L / N) .- A*sin(L*α/N) .- A^3*k^2*sin(L*α/N) .- A^4*k^3 / 3 * sin(2*L*α/N) for α in 1:N];
Y = [(cos(L * α / N )) * A + 1/6*A^4*k^3*cos(L*α/N) .+ (A^2*k / 2).+ A^4*k^3 * 1/2 for α in 1:N];
ϕ = [sqrt(9.81/k) * A  * exp.(k*Y[α]) * sin(k*X[α]) for α in 1:N];




vars = matread(projectdir()*"/data/ClamondICSteepest.mat")
X = vec(vars["X"])
Y = vec(vars["Y"]);
ϕ = vec(vars["F"]);

sol = runSim(N,X,Y,ϕ,Δt,tf,L,smoothing=smoothed,alg = alg)


# Get X, Y, ϕ values

t = 0:Δt:1
xvals = zeros(length(t),N);
yvals = zeros(length(t),N);
ϕvals = zeros(length(t),N);
for i ∈ 1:length(t)
    xvals[i,:] = sol(t[i])[1:N]
    yvals[i,:] = sol(t[i])[N+1:2*N]
    ϕvals[i,:] = sol(t[i])[2*N+1:3*N]
end




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
        xvals = mod.(real.(sol(i)[1:N]),2π)
        yvals = imag.(sol(i)[1:N])
        ϕvals = sol(i)[N+1:2*N]

        scatter(xvals,yvals, legend = false,
         framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1,
          dpi = 300, xlabel=L"x \,(m)",ylabel=L"z \,(m)", title= @sprintf("Time: %.1f s", i),
          xlims=(0,L),ylims = (-1.9,1.9),aspect_ratio=1)
        #plot!([maxVec[Int(t÷Δt + 1)],maxVec[Int(t÷Δt + 1)]],[-2,2],linewidth=3)
    end every interval
    gif(anim, projectdir()*"/plots/RK4noSmooth.gif", fps=fps)
end
visualize(20, 30)

energy, MWL_check = computeEnergy(sol,N,Δt,tf)
plot(0:Δt:tf,energy .- energy[1],xlabel=L"t \, (s)",ylabel =L"E - E_0 \quad (J/kg)",legend=false,title= L"E_0: \quad "*@sprintf("%.4f ", energy[1])*L"(J/kg)")
plot(0:Δt:tf,MWL_check,xlabel=L"t \, (s)",ylabel ="MWL"*L" \quad (m)",legend=false)

# Run simulations for range of steepnesses
for A ∈ Arange
    X = [(α * L / n) - 0*A*sin(k*α/N) - 0*A^3*k^2*sin(k*α/N) - 0*A^4*k^3 / 3 * sin(2*k*α/N) for α in 1:N]
    Y = [(cos(k * α / N )) * A + 0*1/6*A^4*k^3*cos(k*α/N) .+ 0*(A^2*k / 2).+ 0*A^4*k^3 * 1/2 for α in 1:N]
    ϕ = [sqrt(9.81/k) * A   * sin(k*X[α]) for α in 1:N]
    sol = runSim(N,X,Y,ϕ,Δt,tf,smoothing=smoothed)
end