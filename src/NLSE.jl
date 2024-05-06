using DrWatson
@quickactivate "Castawave"
using DifferentialEquations
using LinearAlgebra
using SparseArrays
using Plots


# Purpose of this script is to simulate the NLSE, which can be written as
# i (Fₜ + 1/2 Fₐ) - ϵ(1/8 Fₐₐ + |F|² F) = 0

# First need to discretize our system
function finite_difference_2nd(n::Int, h::Float64,ϵ)
    main_diag = fill(-2.0, n)
    off_diag = fill(1.0, n-1)
    
    matrix = spdiagm(-1=>off_diag, 0=>main_diag,1=>off_diag)

    # Adjust the first and last rows for periodic boundary conditions
    matrix[1, end] = 1.0
    matrix[end, 1] = 1.0
    return h^(-2)*matrix * ϵ^2
end

function finite_difference_1st(n::Int,h::Float64,ϵ)
    top_diag = fill(1.0,n-1)
    bot_diag = fill(-1.0,n-1)

    matrix = spdiagm(-1=>bot_diag,1=>top_diag)
    # periodic BC 
    matrix[1,end] = -1.0
    matrix[end,1] = 1.0
    return matrix / (2 *h) * ϵ
end

function finite_difference_3rd(n::Int,h::Float64,ϵ)
    first_top_diag = fill(-1.0,n-1)
    second_top_diag = fill(1/2,n-2)
    first_bot_diag = fill(1.0,n-1)
    second_bot_diag = fill(-1/2,n-2)

    matrix = spdiagm(-1=>first_bot_diag,1=>first_top_diag,-2=>second_bot_diag,2=>second_top_diag)
    # periodic BC 
    matrix[1,end] = 1.0
    matrix[end,1] = -1.0
    matrix[end-1,1] = 1/2
    matrix[end,2] = 1/2
    matrix[1,end-1]= -1/2
    matrix[2,end] = -1/2
    return matrix / (h)^3 * ϵ^3
end


# Set up function to take derivative 
function NLSE(dF, F, p,t)
    ϵ = p[1]
    D1 = p[2]
    D2 = p[3]
    D3 = p[4]
    #dF .= -im*ϵ*(1/8 * D2 * F .+ abs.(F).^2 .* F) .- 1/2 * D1*F .+ ϵ^2 *(1/8 * D3 * F .- abs.(F).^2 .* D1*F .+ 1/2 * F .* D1 * abs.(F).^2) 
    dF .= -im*ϵ*(1/8 * D2 * F .+ abs.(F).^2 .* F) .- 1/2 * D1*F .+ ϵ^2 *(1/8 * D3 * F .- 1/2* F.* (conj.(F).*D1*F .- F .* D1*conj.(F)) )
end


# Write ODE problem
n = 200
ϵ = 0.1
A = range(-1,stop=20,length=n)
h = A[2]-A[1]
tf = 70.0;

# Initial condition
F_initial = 0.2*sech.((A.-2)).^2 .* exp.(im.* A./ϵ)
plot(A,real.(F_initial))
plot!(A,imag.(F_initial))

params = (ϵ, finite_difference_1st(n,h,ϵ),finite_difference_2nd(n,h,ϵ),finite_difference_3rd(n,h,ϵ))
prob = ODEProblem(NLSE,F_initial,tf,params);

sol = solve(prob,Vern7(),reltol=1e-6)

# Animate solution 
t = 0:0.2:tf
using LaTeXStrings
using Printf
gr()
anim = @animate for i ∈ t
    Fvals = sol(i)

    plot(a,abs.(Fvals), legend = false,
        framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1,
        dpi = 300, xlabel=L"a \,(m)",ylabel=L"F \,(m)", title= @sprintf("T: %.2f s", i),
        xlims=(0,maximum(a)),ylims = (-1.5,1.5))
    plot!(a,real.(Fvals))
    #plot!([maxVec[Int(t÷Δt + 1)],maxVec[Int(t÷Δt + 1)]],[-2,2],linewidth=3)
end every 1
gif(anim, projectdir()*"/plots/NLSE.gif", fps=20)

# Compute drifts
# To second order
U2nd = zeros(n,length(sol.t));
U3rd = zeros(n,length(sol.t));
U4th = zeros(n,length(sol.t));

for i ∈ 1:n
    Fplus1 = sol[mod1(i+1,n),:]
    Fminus1 = sol[mod1(i-1,n),:]
    F = sol[i,:]
    Fₐ = (Fplus1 .- Fminus1)./(2*h)
    # momentum
    momentum = conj(F).*Fₐ .- F.*conj(Fₐ)
    U2nd[i,:] = ϵ^2 * abs.(F).^2
    U3rd[i,:] = -ϵ^3 * im * momentum .* 3/2
    U4th[i,:] = 2* ϵ^4 * abs.(F).^4
end
using Trapz
plot(a[500:1500],trapz(sol.t,U2nd[500:1500,:]))
plot(a[500:1500],trapz(sol.t,U3rd[500:1500,:]))
plot(a[500:1500],trapz(sol.t,U4th[500:1500,:]))
plot(a,abs.(sol(tf)).^2)
plot!(a,abs.(sol(0)).^2)
plot!(a,abs.(sol(35).^2))

plot(sol.t,trapz(1:n,U2nd')) #Conserved
plot(sol.t,trapz(1:n,U3rd')) # Kinda conserved
plot(sol.t,trapz(1:n,U4th')) # Not conserved
plot(sol.t,trapz(1:n,U3rd') .+1/4 *trapz(1:n,U4th'))


plot(abs.(sol.u[1]))
particle = 450
loc = round(a[particle],sigdigits=3)
plot(sol.t,abs.(sol[particle,:]).^2,xlabel="T",ylabel = "F",title = "A = $loc")
using Trapz
trapz(sol.t[1:4000],(abs.(sol[particle,:]).^2)[1:4000])
trapz(sol.t[5000:end],(abs.(sol[particle,:]).^2)[5000:end])
# Integrated drift to second order is constant 

# At third order
# Calculate Fₐ
Fplus1 = sol[particle+1,:]
Fminus1 = sol[particle-1,:]
F = sol[particle,:]
Fₐ = (Fplus1 .- Fminus1)./(2*h)
# momentum
momentum = conj(F).*Fₐ .- F.*conj(Fₐ)
plot(sol.t,abs.(im.*momentum))
# Approximation for |F|_t^2
F² = abs.(F).^2
F²ₜ = diff(F²)./(diff(sol.t))
plot(sol.t[1:end-1],1/2*ϵ^3* F²ₜ.*sol.t[1:end-1].^2)
plot(sol.t,abs.(F).^4)
trapz(sol.t[1:4000],(abs.(F).^4)[1:4000])
trapz(sol.t[5000:end],(abs.(F).^4)[5000:end])
plot(sol.t[1:end-1],F²ₜ .-(abs.(F).^4)[1:end-1])


