using DrWatson
@quickactivate "Castawave"
using LaTeXStrings
using Printf
using Plots
plotlyjs()

# Plot Gerstner (x,y,t) in 3D spacetime
A = 0.3;
k = 1;
g = 1;

# x(a,t) = a .+ A*sin.(k*(a .- (1 - A^2 *k^2)*t)) .+ A^2*k^2 * t
# y(a,t) = -A*cos.(k*(a .- (1 - A^2*k^2)*t))
# metric(a,t) = 1 + A^2*k^2 + A*k*cos.(k*(a .- t))

x(a,t) = a .+ A*sin.(k*(a .- t)) 
y(a,t) = -A*cos.(k*(a .-t))
metric(a,t) = 1 + A^2*k^2 + A*k*cos.(k*(a .- t))

arange = range(-2π,stop=2π,length=100)
trange = range(0,stop=10,length=100)

Plots.surface(x.(arange',trange),repeat(trange,1,length(arange)),y.(arange',trange),surfacecolor=y.(arange',trange),
xlabel = "x(m)",ylabel = "ct (m)",zlabel="y(m)",zlims=(-5,5),colorbar=true,cbar_title="Elevation",aspect_ratio=1)

Plots.plot!(x.(0,trange),trange,y.(0,trange),linewidth=6,linecolor=:red,label="Lagrangian trajectory")

Plots.plot!(x.(arange',trange)',repeat(trange,1,length(arange))',y.(arange',trange)',linewidth=2,legend=false)
Plots.plot!(x.(arange',trange),repeat(trange,1,length(arange)),y.(arange',trange),linewidth=2,legend=false)

Plots.plot(arange,trange,metric.(arange',trange))


# Does the EoM give the right path?
using DifferentialEquations

Xfunc(s,t) = s .+ A*sin.(s .- t)
Yfunc(s,t) = -A*cos.(s .- t)
XS(s,t) = 1 .+ A*cos.(s.-t)
YS(s,t) = A*sin.(s.-t)
Xt(s,t) = -A*cos.(s.-t)
Yt(s,t) = -A*sin.(s.-t)
XSS(s,t) = -A*sin.(s.-t)
YSS(s,t) = A*cos.(s.-t)
XSt(s,t) = A*sin.(s.-t)
YSt(s,t) = -A*cos.(s.-t)
Xtt(s,t) = -A*sin.(s.-t)
Ytt(s,t) = A*cos.(s.-t)

# function f(du,u,p,τ)
#     du[1] = u[2]
#     du[2] = A *(g .+ A^2 .*(u[4] .- u[2]).^2 .+ A *(g .+ (u[4] .- u[2]).^2) .* cos.(u[3] .- u[1])) .* sin.(u[3].-u[1]) ./
#     (1 .+ 2*A^2 .+ 2*A*cos.(u[3].-u[1]) .- A^2 .* cos.(u[3].-u[1]).^2)
#     du[3] = u[4]
#     du[4] =  A .*(-g .+ (1 + A^2).* (u[4] .- u[2]).^2 .+ A *g *cos.(u[3] .- u[1])).* sin.(u[3] .- u[1]) ./
#     (1 .+ 2*A^2 .+ 2*A*cos.(u[3].-u[1]) .- A^2 .* cos.(u[3].-u[1]).^2)
# end

# function f(du,u,p,τ)
#     du[1] = u[2]
#     du[2] = A*(g *(2* u[4] .- 3 .* u[2]).* u[2] .+ 
#     A^2 .* ((u[4] .- u[2]).^2 .+ g *(-1 .+ 2 *u[4].* u[2] .- 3* u[2].^2)) .+ 
#     A *((u[4] .- u[2]).^2 .+ g *(-1 .+ 4 *u[4] .*u[2] .- 6 .*u[2].^2)).* cos.(u[3] .- u[1])) .*sin.(u[3] .- u[1]) ./
#     ((1 .+ 2*A^2 .+ 2*A*(1 + g + A^2 *g).* cos.(u[3].-u[1]) .+ A^2 .*(4*g - 1) .* cos.(u[3].-u[1]).^2))
#     du[3] = u[4]
#     du[4] =  -A *(-A* g *(2 *u[4] .- 3 *u[2]).* u[2] .*(A .+ cos.(u[3] .- u[1])) .+ (g .- (u[4] .- u[2]).^2) .*(1 .+ A^2 .+ 
#     2 *A *g* cos(u[3] .- u[1]))).* sin(u[3] .- u[1]) ./
#     (1 .+ 2*A^2 .+ 2*A*(1 + g + A^2 *g).* cos.(u[3].-u[1]) .+ A^2 .*(4*g - 1) .* cos.(u[3].-u[1]).^2)
# end

function f(du,u,p,t)
    du[1] = u[2]
    du[2] = -(A .* (-1 .+ g .+ 2 * u[2] .- u[2].^2).*sin.(u[1].-t))./(1 .+ A^2 .+ 2*A*cos.(u[1].-t))
end



# initial = [0,1,0,0]
# tauspan = (0,3)
# taurange = range(tauspan[1],stop=tauspan[2],step=0.01)

initial = [0,0]
tspan = (0,trange[end])



trajectory = ODEProblem(f, initial, tspan)
sol = solve(trajectory, Vern7(), reltol = 1e-10);


solution = sol.(trange);
svals = zeros(length(trange));
for i ∈ 1:length(collect(trange))
    svals[i] = solution[i][1]
end

# mean curvature 
meanCurv(s,t) = 2*sqrt(2)*A*(A .- cos.(s .- t))./ (2 .+ 3*A^2 .- A.^2 .* cos.(2 .*(s.-t)).- 4*A*cos.(s.-t).^2)

acc(s,t) = 

plotlyjs()
Plots.surface(x.(arange',trange),repeat(trange,1,length(arange)),y.(arange',trange),surfacecolor=meanCurv.(arange',trange),
xlabel = "x(m)",ylabel = "ct (m)",zlabel="y(m)",colorbar=true,zlims=(-2,2),cbar_title="Mean curvature",aspect_ratio=1)

Plots.plot!(x.(0,trange),trange,y.(0,trange),linewidth=6,linecolor=:red,label="Lagrangian trajectory")
plot!(Xfunc.(svals,trange),trange,Yfunc.(svals,trange),color=:green,linewidth=6,legend=false)



 