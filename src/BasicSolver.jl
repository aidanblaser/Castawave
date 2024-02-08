using DrWatson
@quickactivate "Castawave"
using LaTeXStrings
using Printf
using MAT
using Trapz

include(projectdir()*"/src/MainSolver.jl")

n = 128
A = 0.3
Δt = 0.01
tf = 10.
L = 2π;
k = 1;
h = 0;
smoothing = false;

X = [(α * L / n) - A*sin(k*α*L/n) - A^3*k^2*sin(k*α*L/n) - A^4*k^3 / 3 * sin(2*k*α*L/n) for α in 1:n]
Y = [(cos(k * α*L / n )) * A + 1/6*A^4*k^3*cos(k*α*L/n) .+ (A^2*k / 2) .+ A^4*k^3 * 1/2 for α in 1:n]
ϕ = [sqrt(GRAVITY/k) * A * exp.(k*Y[α]) * sin(k*X[α]) for α in 1:n]
# Use if you want a generic initial condition
# 
# using MAT 
if false
    vars = matread(projectdir()*"/data/ClamondICn128.mat")
    X = vec(vars["X"]);
    Y = vec(vars["Y"]);
    ϕ = vec(vars["F"]);
end

n = length(X)

scatter(X,Y)


MWL = sum(Y)/n
MWL = trapz(X,Y)/(2π)
sol= runSim(n, X, Y, ϕ, Δt, Float64(tf),L,h)



# END OF SIMULATION CODE (remaining code are tests or modified methods)

function visualize(interval::Int, fps::Int)
    anim = @animate for t = 0:Δt:tf
        xvals = sol(t)[1:n]
        yvals = sol(t)[n+1:2*n]
        ϕvals = sol(t)[2*n+1:3*n]

        scatter([xvals],[yvals], label = "Timestepped", legend = :bottomright,
         framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1,
          dpi = 300, xlabel=L"x \,(m)",ylabel=L"z \,(m)", title= @sprintf("Time: %.1f s", t),
          xlims=(0,L),ylims = (-1.9,1.9),aspect_ratio=1)
    end every interval
    gif(anim, projectdir()*"/plots/RK4noSmooth.gif", fps=fps)
end
visualize(1, 30)

function simpsons_rule_periodic(X, Y)
    n = length(X)
    h = diff(X)
    sum = 0   
    for i in 1:2:n-3
        sum += (h[i] / 3) * (Y[i] + 4*Y[i+1] + Y[i+2])
    end
    #do the last point periodically
    sum += (h[end] / 3) * (Y[n-1] + 4*Y[n] + Y[1])
    
    return sum
end


function computeEnergy()
    energy = []
    MWL_check = []
    for t ∈ 0:Δt:tf
        x = sol(t)[1:n]
        y = sol(t)[n+1:2*n]
        ẋ = sol(t,Val{1})[1:n]
        ẏ = sol(t,Val{1})[n+1:2*n]
        ϕ = sol(t)[2*n+1:end]
        R = x .+ im.*y
        Ω,r,θ = conformalMap(R)
        Ω_ξ = DDI1(Ω, n)
        R_ξ = (im ./ Ω) .* Ω_ξ
        q = ẋ .+ im.*ẏ
        ϕ_n = imag.(q .* abs.(R_ξ) ./ R_ξ)
        # Other way of getting KE from Balk 
        B_ξ = ẋ.*imag.(R_ξ) .- ẏ.*real.(R_ξ)
        integrand = -1/2 * ϕ .* B_ξ
        KE = simpsons_rule_periodic(x,ϕ.* ϕ_n/2)
        #KE = trapz(1:n,integrand)
        #println(KE)
        PE = simpsons_rule_periodic(x,GRAVITY/2 * (y).^2)
        append!(MWL_check,simpsons_rule_periodic(x,y)/(2π)) # Eulerian MWL
        #append!(MWL_check,sum(y)/n) # Lagrangian MWL
        append!(energy,KE + PE)
    end
    return energy, MWL_check
end

energy, MWL_check = computeEnergy()
meanE = sum(energy/(tf/Δt))
plot(0:Δt:tf,energy,xlabel=L"t \, (s)",ylabel ="Energy"*L" \quad (J/kg)",legend=false)
plot(0:Δt:tf,MWL_check,xlabel=L"t \, (s)",ylabel ="MWL"*L" \quad (m)",legend=false)
savefig(projectdir()*"/plots/Energy")