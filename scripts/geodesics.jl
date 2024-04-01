function spaceCurvature(xvals,yvals,t,n)
    curvature_space = copy(xvals)
    Ωξ = zeros(Complex,(length(t),n))
    Ωξξ = zeros(Complex,(length(t),n))
    Rξ = zeros(Complex,(length(t),n))
    Rξξ = zeros(Complex,(length(t),n))
    Xξ = copy(xvals)
    Yξ = copy(xvals)
    Xξξ = copy(xvals)
    Yξξ = copy(xvals)
    R = xvals .+ im.* yvals
    Ω = exp.(- im * R);
    for i = 1:length(t)
        Ωξ[i,:] = DDI1(Ω[i,:],n,0)
        Ωξξ[i,:]= DDI2(Ω[i,:],n,0)
        Rξ[i,:] = (im ./ Ω[i,:]) .* Ωξ[i,:]
        Rξξ[i,:] = im.*(Ω[i,:].*Ωξξ[i,:] .- Ωξ[i,:].^2)./(Ω[i,:].^2)
        Xξ[i,:] =  real.(Rξ[i,:])
        Yξ[i,:] =  imag.(Rξ[i,:])
        Xξξ[i,:] =  real.(Rξξ[i,:])
        Yξξ[i,:] =  imag.(Rξξ[i,:])
        curvature_space[i,:] = (Xξ[i,:].*Yξξ[i,:] .- Yξ[i,:].*Xξξ[i,:])./((Xξ[i,:].^2 .+ Yξ[i,:].^2).^(3/2))
    end
    return curvature_space
end

using CubicSplines
function timeCurvature(xvals,yvals,t,n)
    curvature_time = zeros(length(t),n)
    Ωξ = zeros(Complex,(length(t),n))
    Ωξξ = zeros(Complex,(length(t),n))
    Rξ = zeros(Complex,(length(t),n))
    Rξξ = zeros(Complex,(length(t),n))
    Xξ = zeros(length(t),n)
    Yξ = zeros(length(t),n)
    Xξξ = zeros(length(t),n)
    Yξξ = zeros(length(t),n)
    Xξt = zeros(length(t),n)
    Yξt = zeros(length(t),n)
    Xt = zeros(length(t),n)
    Yt = zeros(length(t),n)
    Xtt = zeros(length(t),n)
    Ytt = zeros(length(t),n)
    R = xvals .+ im.* yvals
    Ω = exp.(- im * R);
    for i = 2:length(t)-1
        Ωξ[i,:] = DDI1(Ω[i,:],n,0)
        Ωξξ[i,:]= DDI2(Ω[i,:],n,0)
        Rξ[i,:] = (im ./ Ω[i,:]) .* Ωξ[i,:]
        Rξξ[i,:] = im.*(Ω[i,:].*Ωξξ[i,:] .- Ωξ[i,:].^2)./(Ω[i,:].^2)
        Xξ[i,:] =  real.(Rξ[i,:])
        Yξ[i,:] =  imag.(Rξ[i,:])
        Xξξ[i,:] =  real.(Rξξ[i,:])
        Yξξ[i,:] =  imag.(Rξξ[i,:])
        Xt[i,:] = (xvals[i+1,:].-xvals[i,:])./(t[i+1]-t[i])
        Yt[i,:] = (yvals[i+1,:].-yvals[i,:]./(t[i+1]-t[i]))
        Xtt[i,:] = (xvals[i+1,:] .- 2*xvals[i,:].+xvals[i-1,:])./((t[i+1]-t[i])^(3/2))
        Ytt[i,:] = (yvals[i+1,:] .- 2*yvals[i,:] .+ yvals[i-1,:])./((t[i+1]-t[i])^(3/2))
        Xξt[i,:] = (Xξ[i+1,:].-Xξ[i,:])./(t[i+1]-t[i])
        Yξt[i,:] = (Yξ[i+1,:].-Yξ[i,:])./(t[i+1]-t[i])
    end
    curvature_time = (Ytt .- 2*(Yξt.*Xt)./Xξ .+ Yξξ.*Xt.^2 ./(Xξ.^2) .+ Yξ.*(2*(Xξt.*Xt)./(Xξ.^2) .+ Xξξ.*Xt.^2 ./(Xξ.^2) .- Xtt./Xξ))./(((Yt .- Yξ.*Xt./Xξ).^2 .+ 1).^(3/2))
    return curvature_time, Xξ,Yξ,Xt,Yt,Xξt,Yξt,Xtt,Ytt,Xξξ,Yξξ
end

function Scurvature(xvals,yvals,t,n)
    curvature_s = zeros(length(t),n)
    Xt = zeros(length(t),n)
    Yt = zeros(length(t),n)
    Xtt = zeros(length(t),n)
    Ytt = zeros(length(t),n)
    for i = 2:length(t)-1
        Xt[i,:] = (xvals[i+1,:].-xvals[i,:])./(t[i+1]-t[i])
        Yt[i,:] = (yvals[i+1,:].-yvals[i,:]./(t[i+1]-t[i]))
        Xtt[i,:] = (xvals[i+1,:] .- 2*xvals[i,:].+xvals[i-1,:])./((t[i+1]-t[i])^(2))
        Ytt[i,:] = (yvals[i+1,:] .- 2*yvals[i,:] .+ yvals[i-1,:])./((t[i+1]-t[i]^(2)))
    end
    curvature_s = sqrt.(((1 .+ Xt.^2 .+Yt.^2).*Xtt .- Xt.^2 .* Xtt .- Xt .* Yt .* Ytt).^2 .+
    (Xt.*Xtt .+ Yt.*Ytt).^2 .+ ((1 .+ Xt.^2 .+ Yt.^2).* Ytt .- Yt.^2 .* Ytt .- Yt.*Xt.*Xtt).^2)./((Xt.^2 .+ Yt.^2 .+ 1).^(3/2))
    return curvature_s
end


time_curvature,Xξ,Yξ,Xt,Yt,Xξt,Yξt,Xtt,Ytt,Xξξ,Yξξ = timeCurvature(xvals,yvals,0:Δt:tf-offset,n)
s_curvature = Scurvature(xvals,yvals,0:Δt:tf-offset,n)
space_curvature = spaceCurvature(xvals,yvals,0:Δt:tf-offset,n)

curvature_capped = min.(max.(space_curvature, -10), 10)
time_curv_capped = min.(max.(time_curvature,-0.1),0.1)
scurv_capped = min.(max.(s_curvature,-0.1),0.5)

surf = surface(xvals,repeat(0:Δt:tf-offset,1,n),yvals,surfacecolor = curvature_capped,
xlabel="x (m)",ylabel="t (s)",zlabel="y (m)",zlims=(-2,2),cbar_title="Spatial Curvature")
plot!(xvals[:,120],0:Δt:tf-offset,yvals[:,120],linewidth=5,linecolor=:red,legend=false)

surf = surface(xvals,repeat(0:Δt:tf-offset,1,n),yvals,surfacecolor = time_curv_capped,
xlabel="x (m)",ylabel="t (s)",zlabel="y (m)",zlims=(-2,2),cbar_title="Time Curvature")
plot!(xvals[:,120],0:Δt:tf-offset,yvals[:,120],linewidth=5,linecolor=:red,legend=false)

surf = surface(xvals,repeat(0:Δt:tf-offset,1,n),yvals,surfacecolor = scurv_capped,
xlabel="x (m)",ylabel="t (s)",zlabel="y (m)",zlims=(-2,2),cbar_title="S Curvature")
plot!(xvals[:,120],0:Δt:tf-offset,yvals[:,120],linewidth=5,linecolor=:red,legend=false)


# Plot values
plotlyjs()
time = 0:Δt:tf-offset
particle = 120
plot(time,space_curvature[:,particle],legend=:bottomleft, label = "Spatial Curvature",ylims=(-100,100))
plot!(time, (Xξ[:,particle].*Xtt[:,particle] .+ Yξ[:,particle].*(Ytt[:,particle] .+ 9.81))./(Xξ[:,particle].^2 .+ Yξ[:,particle].^2)/30, label= "Tangential Acceleration")
#plot!(time,(Yξ[:,particle].*(Ytt[:,particle] .+ 9.81))./(Xξ[:,particle].^2 .+ Yξ[:,particle].^2)/30)
plot!(time,(Xξ[:,particle].*Xtt[:,particle]./(Xξ[:,particle].^2 .+ Yξ[:,particle].^2)),label = "Horizontal Component of Tangential Acceleration")
plot!(time,Xtt[:,particle].*10,label="Horizontal Acceleration")
plot!(time,(Xξ[:,particle].*Xt[:,particle] .+ Yξ[:,particle].*Yt[:,particle]).*10,label="Kappa term")

plot(time,Xtt[:,particle],label="Horizontal Acceleration",legend=:topleft)
plot!(time, -9.81*Yξ[:,particle].*yvals[:,particle].*25,label="Yξ * Y")
plot!(time, -9.81*Yξ[:,particle]./Yt[:,particle])
plot!(time,Xξ[:,particle]*10,label ="Xξ")

plot(time,(Xtt[:,particle].*Yξ[:,particle].- Ytt[:,particle].*Xξ[:,particle].-9.81*Xξ[:,particle])./(Xξ[:,particle].^2 .+ Yξ[:,particle].^2),legend=:topleft,label="pb / J")
plot!(time,-1.9*9.81*maximum(yvals,dims=2).*Xξ[:,particle]./(Xξ[:,particle].^2 .+ Yξ[:,particle].^2),label="-g Xξ * waveheight / (Xξ^2 + Yξ^2)")
plot(time,(Xtt[:,particle].*Yξ[:,particle])./(Xξ[:,particle].^2 .+ Yξ[:,particle].^2))
plot!(time,-1.9*9.81*maximum(yvals,dims=2).*Xξ[:,particle]./ Yξ[:,particle ]./(Xξ[:,particle].^2 .+ Yξ[:,particle].^2))
plot(time,Yξ[:,particle])

# Normal component of accelerations
plot(time,(-Xtt[:,particle].*Yξ[:,particle].+ Ytt[:,particle].*Xξ[:,particle])./sqrt.(Xξ[:,particle].^2 .+ Yξ[:,particle].^2),legend=:topleft,label="Normal Acceleration")
plot!(time, - Xξ[:,particle]./sqrt.(Xξ[:,particle].^2 .+ Yξ[:,particle].^2) ,label="Normal Component of Gravity (g's)")


# Tangential component of accelerations
plot(time, (Xtt[:,particle].*Xξ[:,particle].+ Ytt[:,particle].*Yξ[:,particle])./(sqrt.(Xξ[:,particle].^2 .+ Yξ[:,particle].^2)),legend=:topleft,label="Tangential Acceleration")
plot!(time, -Yξ[:,particle] ./ sqrt.(Xξ[:,particle].^2 .+ Yξ[:,particle].^2),label="Tangential component of gravity (g's)")
plot(time,(Xtt[:,particle].*Xξ[:,particle].+ Ytt[:,particle].*Yξ[:,particle] .+9.81* Yξ[:,particle]))

plot(time,Xtt[:,particle],legend=:topleft)
plot!(time,-9.81*100*maximum(yvals,dims=2).*Xξ[:,particle])