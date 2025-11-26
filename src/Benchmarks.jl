using DrWatson
@quickactivate "Castawave"
using LaTeXStrings
using Printf
using MAT
using Roots

include(projectdir()*"/src/MainSolver.jl")

# Initial parameters
n = 128
A = 0.3
Δt = 0.01
tf = 9
L = 2π;
k = 1;
h = 0.0;

# Build Stokes Wave initial condition
X = [(α * L / n) - A*sin(k*α*L/n) - A^3*k^2*sin(k*α*L/n) - A^4*k^3 / 3 * sin(2*k*α*L/n) for α in 1:n]
Y = [(cos(k * α*L / n )) * A + 1/6*A^4*k^3*cos(k*α*L/n) .+ (A^2*k / 2).+ A^4*k^3 * 1/2 for α in 1:n]
ϕ = [sqrt(GRAVITY/k) * A * exp.(k*Y[α]) * sin(k*X[α]) for α in 1:n]

if true
    vars = matread(projectdir()*"/data/ClamondICn128.mat")
    X = vec(vars["X"]);
    Y = vec(vars["Y"]);
    ϕ = vec(vars["F"]);
end
n = length(X)

# p = scatter(X,Y,aspect_ratio=1,xlabel="x (m)",ylabel="y (m)",title="Initial condition")
# savefig(p,projectdir()*"/plots/IC")

# Analytic solution at surface
U(A,k,g=GRAVITY) = (A^2*k^2 + A^4*k^4/2 + A^6*k^6*37/24 + A^8*k^8*1739/360)*sqrt(g/k)*(1 + A^2*k^2/2 + A^4*k^4/8 + A^6*k^6/16)
θ(a,t,A,k,g=GRAVITY) = k*(a .- (sqrt(g/k).*(1 + 1/2*A^2*k^2 + 1/8*A^4*k^4 .+ A^6*k^6/16) .- U(A,k)).*t)
x(a,t,A,k,g=GRAVITY) = a .+ U(A,k).*t .- A*sin.(θ.(a,t,A,k)) .- A^3*k^2*sin.(θ.(a,t,A,k)) .- 
 1/3*A^4*k^3*sin.(2*θ.(a,t,A,k)) .- 47/24*A^5*k^4*sin.(θ.(a,t,A,k)) .- 5/72*A^5*k^4*sin.(3*θ(a,t,A,k)) .-
 17/18*A^6*k^5*sin.(2*θ(a,t,A,k)) .- 1/80*A^6*k^5*sin.(4*θ(a,t,A,k))
y(a,t,A,k,g=GRAVITY) = 1/2*A^2*k .+ A^4*k^3/2 .+23/24*A^6*k^6 .+ 
A*cos.(θ.(a,t,A,k)) .+ A^4*k^3*cos.(2*θ.(a,t,A,k)) .- 1/24*A^5*k^4*cos.(θ.(a,t,A,k)).+ 
11/36*A^6*k^5*cos.(θ.(a,t,A,k)) .+ 1/120*A^6*k^5*cos.(4*θ.(a,t,A,k))

arange = range(0,step=2π/(10*n),length=10*n)

# Run simulation
Xfull, Yfull, ϕfull, t= @time runSim(n, X, Y, ϕ, Δt, Float64(tf),L,h)

plot(Xfull[end,:],Yfull[end,:])


anim = @animate for i ∈ 1:length(t)
    plot(x.(arange,t[i],A,k),y.(arange,t[i],A,k),label="Analytic",linewidth=3,legend = :bottomright,
    framestyle= :box,
     dpi = 300, xlabel=L"x \,(m)",ylabel=L"z \,(m)", title= @sprintf("Time: %.1f s", t[i]),
     xlims=(0,L),ylims = (-1.9,1.9),aspect_ratio=1)
    scatter!(Xfull[i,:],Yfull[i,:], label = "Simulation",markerstrokewidth=0, markersize=2)
end
gif(anim, projectdir()*"/plots/AnalyticComparison.gif", fps=30)




# compute RMSE
rmse = []
for t = 0:Δt:tf
    xvals = sol(t)[1:n]
    yvals = sol(t)[n+1:2*n]
    ϕvals = sol(t)[2*n+1:3*n]

    # find the a-values that correspond to the x-values of numerical output
    avals = []
    for i=1:n
        xdev(a) = x.(a,t,A,k) .- xvals[i]
        append!(avals,find_zero(xdev,xvals[i]))
    end



    append!(rmse,sqrt(sum((yvals .- y(avals,t,A,k)).^2)/n))
end

p = plot(0:Δt:tf,rmse,xlabel="time (s)",ylabel = "RMSE (m)",label = L"\sqrt{\langle (y - y_p)^2\rangle}",legendfont=14)
savefig(p,projectdir()*"/plots/RMSE")

energy, MWL_check = computeEnergy()
plot(0:Δt:tf,energy .- energy[1],xlabel=L"t \, (s)",ylabel =L"E - E_0 \quad (J/kg)",legend=false,title= L"E_0: \quad "*@sprintf("%.4f ", energy[1])*L"(J/kg)")
plot(0:Δt:tf,MWL_check,xlabel=L"t \, (s)",ylabel ="MWL"*L" \quad (m)",legend=false)
savefig(projectdir()*"/plots/MWL")