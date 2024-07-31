using DrWatson
@quickactivate "Castawave"
using LaTeXStrings
using Printf
using MAT
using Trapz
using Plots
plotlyjs()

include(projectdir()*"/src/MainSolver.jl")
include(projectdir()*"/src/ClamondIC.jl")
include(projectdir()*"/src/DoldMono.jl")

N = 256
A = 0.4
Δt = 0.001
tf = 5
L = 2π;
k = 1;
h = 0;
smoothing = false;
alg = Vern9()

X = [(α * L / n) - 0*A*sin(k*α*L/n) - 0*A^3*k^2*sin(k*α*L/n) - 0*A^4*k^3 / 3 * sin(2*k*α*L/n) for α in n÷2:n+n÷2-1]
Y = [(cos(k * α*L / n )) * A + 0*1/6*A^4*k^3*cos(k*α*L/n) .+ 0*(A^2*k / 2).+ 0*A^4*k^3 * 1/2 for α in n÷2:n+n÷2-1]
ϕ = [sqrt(9.81/k) * A   * sin(k*X[α]) for α in 1:n]
# Use if you want a generic initial condition
scatter(X,Y)
scatter!(X,ϕ)
# 
# using MAT 
if true
    vars = matread(projectdir()*"/data/ClamondICn128.mat")
    X = vec(vars["X"]);
    Y = vec(vars["Y"]);
    ϕ = vec(vars["F"]);
end

n = length(X)

scatter(X,Y)


MWL = sum(Y)/n
MWL = simpsons_rule_periodic(X,Y)/(2π)

N = 256;
A = 0.4
X,Y,ϕ,c = getIC(Inf,A,N÷2)
tf = 1;
Δt = 1e-3
L = 2π
h = 0.0
tol = 1e-6
smoothed = false

ZC = X .+ im.*Y
Zshift = ZC.* exp.(im*π/2)
scatter(real(ZC),imag(ZC))

scatter(X,Y)

Xfull, Yfull, ϕfull, t = @time runSim(N, X, Y, ϕ, Δt, Float64(tf),L,h,ϵ = tol,smoothing=smoothed)
t[end]
plot(t[6:end-1],diff(t)[6:end],xlabel="t (s)",ylabel="dt (s)")

sum(Yfull[1,:] .*DDI1(Xfull[1,:],N,2π,1))/N
sum(Yfull[end,:].*DDI1(Xfull[end,:],N,2π,1))/N

# Check relations 
half = length(t)÷2 + 100 
ϕhalf = ϕfull[half,:]
Xhalf = Xfull[half,:]
Yhalf = Yfull[half,:]
xξ = DDI1(Xhalf,N,L,1)
yξ = DDI1(Yhalf,N,0,1)
(ẋ, ẏ, _, _, _, _, _, _, _) = fixedTimeOperations(N, Xhalf, Yhalf, ϕhalf, L, 0.0);
ψξ = ẋ.*yξ .- ẏ.*xξ
ϕξ = DDI1(ϕhalf,N,0,1)
ϕξL = ẋ.*xξ .+ ẏ.*yξ

X = Xhalf
Y = Yhalf
Ω = conformalMap(X .+ im*Y)
R_ξ = DDI1(X,N,L,1) .+ im*DDI1(Y,N,0,1)
R_ξξ = DDI2(X,N,L,1) .+ im*DDI2(Y,N,0,1)
Ω_ξ = DDI1(Ω,N,0,0)
Ω_ξξ = DDI2(Ω,N,0,0)

H = conformalDepth(h)

# The matrix method described in Dold is used to find the normal derivative of the potential.
A, B, ℵ = ABMatrices(Ω, Ω_ξ, Ω_ξξ, N, H)
ϕ_ξ, ϕ_ν = NormalInversion(ϕhalf, A, ℵ, N)
plot(ϕ_ν)
plot!(-ψξ)
plot(ϕ_ν)
plot(ϕ_ν .+ ψξ,label="ϕν + ψξ")
plot!(ϕξL)
plot(ϕ_ξ.-ϕξ)



plot(ϕξ.-ϕξL)
plot(ψξ,label="ψξ")
plot!(ϕξ,label="ϕξ")
plot!(-imag.(hilbert(ϕξ)),label="H(ϕξ)")

# END OF SIMULATION CODE (remaining code are tests or modified methods)
gr()
function visualize(interval::Int, fps::Int)
    anim = @animate for i ∈ enumerate(t)
        scatter(mod.(Xfull[i[1],:],L),Yfull[i[1],:], legend = false,
         framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1,
          dpi = 300, xlabel=L"x \,(m)",ylabel=L"z \,(m)", title= @sprintf("Time: %.1f s", i[2]),
          xlims=(0,L),ylims = (-1.9,1.9),aspect_ratio=1)
    end every interval
    gif(anim, projectdir()*"/plots/RK4noSmooth.gif", fps=fps)
end
visualize(10, 10)

KE,PE,MWL,phasespd,momentum = computeEnergyDold(Xfull,Yfull,ϕfull,t,N,c)

xD,yD,ϕD,tD = monoSim(X,Y,ϕ,10)
tD[end]
sum(yD[1,:] .*DDI1(xD[1,:],N,2π,1))/N
sum(yD[end,:].*DDI1(xD[end,:],N,2π,1))/N

KED,PED,MWLD,phasespdD,momentumD = computeEnergyDold(xD,yD,ϕD,tD,N,c)
plot(tD,momentumD)
plot(tD,MWLD)
plot(tD,KED.+PED)
plot(tD,phasespdD)


# t = [120]
#     x = Xfull[t[1],:]
#     y = Yfull[t[1],:]
#     ϕ = ϕfull[t[1],:]
#     xξ = DDI1(x,N,L,1)
#     yξ = DDI1(y,N,0,1)
#     (ẋ, ẏ, _, _, _, _, _, _, _) = fixedTimeOperations(N, x, y, ϕ, L, 0.0)
#     B_ξ = ẋ.*yξ .- ẏ.*xξ
#     ϕ_ξ = DDI1(ϕ,N,0,1)

# plot(B_ξ,label="B_ξ")
# plot!(ϕ_ξ,label="ϕ_ξ")
# plot!(imag.(hilbert(B_ξ)),label="H(B_ξ)")
# plot(abs.(ϕ_ξ .- imag.(hilbert(B_ξ))))
# abs.(ϕ_ξ .- imag.(hilbert(B_ξ)))


MWL = sum(y.*xξ)/ N
B_ξ = ẋ.*yξ .- ẏ.*xξ

# Compare MWL
gr()
title = "N=$N, Ak=$A, ϵ = 1e-6"
plot(t,abs.(MWL.-MWL[1].-1e-16),xlabel="t (s)",ylabel = "MWL (m)",label="Castawave",yscale=:log10,title=title)
plot!(tD,abs.(MWLD.-MWLD[1].-1e-16),label="Dold")
# Compare energies
plot(t,abs.(KE.+PE .- (KE[1]+PE[1]).-1e-16),xlabel="t (s)",ylabel = "Energy (J/kg)",label="Castawave",yscale=:log10,title=title)
plot!(tD,abs.(KED.+PED.-(KED[1]+PED[1]).-1e-16),label="Dold")
# Compare phasespd
plot(t,abs.(phasespd.-phasespd[1].-1e-16),xlabel="t (s)",ylabel = "c (m/s)",label="Castawave",yscale=:log10,title=title)
plot!(tD,abs.(phasespdD.-phasespdD[1].-1e-16),label="Dold")

gr()
function visualize(interval::Int, fps::Int)
    anim = @animate for i ∈ enumerate(tD)
        scatter(mod.(xD[i[1],:],2π),yD[i[1],:], legend = false,
         framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1,
          dpi = 300, xlabel=L"x \,(m)",ylabel=L"z \,(m)", title= @sprintf("Time: %.1f s", i[2]),
          xlims=(0,L),ylims = (-1.9,1.9),aspect_ratio=1)
    end every interval
    gif(anim, projectdir()*"/plots/RK4noSmooth.gif", fps=fps)
end
visualize(1, 10)


energy, MWL_check = computeEnergy(sol,n,Δt,tf)
plot(0:Δt:tf,energy .- energy[1],xlabel=L"t \, (s)",ylabel =L"E - E_0 \quad (J/kg)",legend=false,title= L"E_0: \quad "*@sprintf("%.4f ", energy[1])*L"(J/kg)")
plot(0:Δt:tf,MWL_check,xlabel=L"t \, (s)",ylabel ="MWL"*L" \quad (m)",legend=false)
savefig(projectdir()*"/plots/Energy")

offset = 0
xvals = zeros(length(0:Δt:tf-offset),n);
yvals = zeros(length(0:Δt:tf-offset),n);
ϕvals = zeros(length(0:Δt:tf-offset),n);

for t = 1:(length(0:Δt:tf-offset))
    xvals[t,:] = sol(Δt*(t-1))[1:n]
    yvals[t,:] = sol(Δt*(t-1))[n+1:2*n]
    ϕvals[t,:] = sol(Δt*(t-1))[2*n+1:3*n]
end

plot(xvals)

plotlyjs()
surface(xvals,repeat(0:Δt:tf-offset,1,n),yvals,
xlabel="x (m)",ylabel="t (s)",zlabel="y (m)",zlims=(-2,2),cbar=false)
plot!(xvals[:,60],0:Δt:tf-offset,yvals[:,60],linewidth=5,linecolor=:red,label="trajectory")
