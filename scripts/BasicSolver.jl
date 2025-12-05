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

N = 256
n = N
A = 0.45
Δt = 0.01
T = 3
L = 2π;
k = 1;
h = 0.0;
smoothing = false;
p = SimulationParameters(L,h,Δt,T,smoothing=true)
X = [(α * L / n) - A*sin(k*α*L/n) for α in 1:n]
Y = [(cos(k * α*L / n )) * A for α in 1:n]
ϕ = [sqrt(9.81/k) * A * sin(k*α*L/n) for α in 1:n]

offset = 100
κ = (DDI2(Yfull[end-offset,:],0).*DDI1(Xfull[end-offset,:],2π) .- DDI2(Xfull[end-offset,:],2π).*DDI1(Yfull[end-offset,:],0))./(DDI1(Xfull[end-offset,:],2π).^2 .+ DDI1(Yfull[end-offset,:],2π).^2).^(3/2)
plot(κ)

Xfull, Yfull, ϕfull, t = @time runSim(X, Y, ϕ,p)
xD,yD,ϕD,tD,uD,vD,KED,PED = @time DoldSim(X, Y, ϕ,T,L)
plotlyjs()
fig = plot(Xfull[end-offset,:],Yfull[end-offset,:],dpi=300,xlims=(4,6),ylims=(0,0.6),xlabel="x (m)",ylabel="y (m)",label="Castawave",title = "A = 0.45 m")
plot!(xD[end-1,:],yD[end-1,:],label="Dold")

t[end-100]
tD[end-2]


plot(xD[1,:],yD[1,:])
plot!(Xfull[1,:],Yfull[1,:])
yD[1,:] .- Yfull[1,:]
#savefig(fig,projectdir()*"/plots/breaking_comparison")
t
Xfull
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

N = 1024;
A = 0.4
X,Y,ϕ,c = getIC(Inf,A,N÷2)
T = 5.0;
Δt = 1e-2
L = 2π
h = 0.0
tol = 1e-6
smoothed = true

ZC = X .+ im.*Y
Zshift = ZC.* exp.(im*π/2)
scatter(real(ZC),imag(ZC))

scatter(X,Y)
scatter!(X,ϕ)

p = SimulationParameters(L,h,Δt,T)

using StatProfilerHTML
StatProfilerHTML.@profilehtml runSim(X, Y, ϕ,p)

plotlyjs()
gr()
plt = plot(Xfull[end,:],Yfull[end,:],legend=false,ylims=(-1,1),xlabel="x (m)",ylabel="y (m)",dpi=300,title="A = 0.45 m Breaker")
#savefig(plt,projectdir()*"/plots/Castawavebreaker")

t
tdesired = collect(range(0,t[end],length=400))
indices = Int64.(ones(length(tdesired)))
for i ∈ 1:length(tdesired)
    indices[i] = argmin(abs.(t .- tdesired[i]))
end


plot(t)


gr()
anim = @animate for i ∈ indices
    time = @sprintf("%0.1f",t[i]);
    scatter(mod.(Xfull[i,:] .+ π,2π) , Yfull[i,:],xlabel=L"x \, \, [m]",ylabel=L"y\, \, [m]",
    title = "t = $time s", titlefont=font(14,"Computer Modern"),
    xlims=(0,2π),ylims=(-1.7,2),guidefont=font(12,"Computer Modern"),tickfont=font(9,"Computer Modern"),
    framestyle= :box,label=false,markersize=2,markerstrokewidth=0,color=:dodgerblue,aspect_ratio=1,dpi=200)
end 
gif(anim, projectdir()*"/plots/breakingVisual.gif", fps = 30)



sum(Yfull[1,:] .*DDI1(Xfull[1,:],N,2π,1))/N
sum(Yfull[end,:].*DDI1(Xfull[end,:],N,2π,1))/N

@load projectdir()*"/data/focusing.jld2" x_o y_o p_o to
Xfull = x_o;
Yfull = y_o;
ϕfull = p_o;
t = to[:,1];
N = 512

# Check relations 
half = length(t)÷2  - 100
half = 1
ϕhalf = ϕfull[half,:]
Xhalf = Xfull[half,:]
Yhalf = Yfull[half,:]
xξ = DDI1(Xhalf,N,L,1)
yξ = DDI1(Yhalf,N,0,1)
(ẋ, ẏ, _, _, _, _, _, _, _) = fixedTimeOperations(N, Xhalf, Yhalf, ϕhalf, L, 0.0);
ψξ = ẋ.*yξ .- ẏ.*xξ
ϕξ = DDI1(ϕhalf,N,0,1)
ϕξL = ẋ.*xξ .+ ẏ.*yξ
plot(ϕξ,ylabel="ϕ_ξ")
sum(ϕξ)


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
plot(ϕ_ν .+ ψξ,label="ϕν + ψξ",xlabel="ξ",ylabel = "ϕν + ẏ xξ - ẋ yξ",title="Deviation from theory")
plot!(ϕξL)
plot(ϕ_ξ.-ϕξ)

ds = xξ.^2 .+ yξ.^2
plot(ds)
dξdt = (imag.(hilbert(ψξ)).-ϕξ)./ ds
plot(dξdt)
plot(ϕξ,label="ϕξ")
plot!(ψξ,label="ψξ")
plot!(-imag.(hilbert(ϕξ)),label="-H(ϕξ)")
plot!(imag.(hilbert(ψξ)),label="H(ψξ)")
plot(imag.(hilbert(ψξ))./ϕξ .* ds)

scatter(X,Y)
plot(imag.(hilbert(ϕ_ξ)).+ψξ)
plot(ẋ)
phix = (ϕξ.*xξ .+ ψξ.*yξ)./ds
plot!(phix)

plot(ψξ)
plot!(imag.(hilbert(ϕ_ξ)))
plot(ϕ_ν,label="ϕ_ν")
plot!(imag.(hilbert(ϕ_ξ )),label="H(ϕξ)")
plot(ϕ_ν .- imag.(hilbert(ϕ_ξ)))
mean(ds)
dsmean = ds./mean(ds);
plot(dsmean)
plot(imag.(hilbert(ϕ_ν)).+ϕξ)

plot(ϕξ.-ϕξL)
plot(ψξ,label="ψξ")
plot!(ϕξ,label="ϕξ")
plot!(ϕ_ν,label="ϕν")
plot!(ψξ .- sqrt.(ds)./(2*c))
plot!(-imag.(hilbert(ϕξ)).- (ds),label="H(ϕξ)")
plot!(imag.(hilbert(ψξ)),label="H(ψξ)")
plot(ψξ .+ imag.(hilbert(ϕξ)) )
plot!()
plot(ϕξ .- imag.(hilbert(ψξ)))
plot!(imag.(hilbert(dsmean)))
plot(xξ)
plot!(sqrt.(xξ.^2 + yξ.^2))
plot(xξ)
plot!(xξ.*sqrt.(xξ.^2 .+ yξ.^2))
plot!(yξ./xξ)
plot(ψξ)
plot!(-imag.(hilbert(ϕξ)) - imag.(hilbert(ϕξ)).*sqrt.(ds))
plot!(-imag.(hilbert(ϕξ)))
plot((ψξ .+ imag.(hilbert(ϕξ))),label="difference")
plot!(ψξ,label="ψξ")
plot!(-imag.(hilbert(ϕξ)),label="-H(ϕξ)")

diff1 = (ψξ .+ imag.(hilbert(ϕξ)));
diff2 = ϕξ .- imag.(hilbert(ψξ));
plot(diff1)
plot!(diff2)
plot(imag.(hilbert(diff1)).+diff2)
plot(diff1.+diff2)

# What can I add to ϕ to make hilbert transform work
plot(ψξ .+ imag.(hilbert(ϕξ .+ ds.*dξdt)))
plot(dξdt,label="dξ/dt")
plot!(X.*N./2π,Y.*10)


plot(ϕξ,label="ϕξ")
plot!(-imag.(hilbert(ϕ_ν)),label="H(ϕν)")
plot!(ϕξ + imag.(hilbert(ϕ_ν)),label="difference")
ψξξ = DDI1(ψξ,N,0,1);
b = A * ψξ .- ψξξ;
ψν = ℵ \ b
plot(ψν)
plot!(ϕξ)
plot(ϕξ.-ψν)
plot(imag.(hilbert(ψν)))
plot!(-ψξ)
plot(imag.(hilbert(ψξ)).-ϕξ)
plot(imag.(hilbert(ϕξ).-ψξ))


# Look at curvature
yξξ = DDI2(Y,N,0,1)
xξξ = DDI2(X,N,2π,1)
κ = (yξξ.*xξ .- xξξ.*yξ)./(ds.^(3/2))
plot(κ)
plot(κ.*ds)


mean(xξ)
scatter(X,Y)

using Trapz
# Try and integrate them 
ψ = zeros(N)
ϕ = zeros(N)
for i ∈ 2:N
    ϕ[i] = trapz((1:i).*2π/N,ϕξ[1:i])
    ψ[i] = trapz((1:i).*2π/N,ψξ[1:i])
end
ϕ = ϕ .- mean(ϕ)
ψ = ψ.-mean(ψ)

plot(ϕ)
plot!(ψ)
plot(ϕ .- imag.(hilbert(ψ)))
plot!(diff1)
plot!(diff2)
plot(sqrt.(ds))


#
function hilbert_circ(ϕ,N)
    xirange = (1:N).*(2π/N)
    val = copy(ϕ)
    for i ∈ 1:N
        integrand = ϕ[i]*cot.((ϕ.-ϕ[i])./2)
        integrand[i] = 0
        val[i] = 1/(2π)*trapz(xirange,integrand)
    end
return val
end







# Energy relationship 
KE = sum(0.5* ϕhalf.*ϕ_ν)
A = 0.3
KE_theory = c/2 * (A^2 / 2 - A^4 / 4 - 35* A^6 / 48)*sqrt(9.81)*2π
PE = sum(0.5*Yhalf.^2 .*xξ)
PE_theory = 2π.*(A^2/4 - A^4/8 - 19* A^6/6)

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
using JLD2
jldsave(projectdir()*"/data/shallowOverturnN1024h06A02.jld2",x = Xfull,y= Yfull,ϕ= ϕfull,t= t,c= c,h= h,N= N,A=A) 

using HDF5
fid=h5open("shallowOverturnN1024h09A03.h5", "w")
fid["X"]= Xfull;
fid["Y"]=Yfull;
fid["Phi"]=ϕfull;
fid["N"]=N;
fid["t"]=t;
fid["A"]=A;
fid["h"] = h;
fid["c"] = c;
close(fid)

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

