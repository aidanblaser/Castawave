using DrWatson
@quickactivate "Castawave"
using LaTeXStrings
using Printf
using MAT
using BenchmarkTools
using JLD2

#=
This script is for testing certain properties of waves. 
Those include:
    - Spectra over time
    - Runtime
    - Energy changes
    - MWL changes
=#

include(projectdir()*"/src/MainSolver.jl")


# Set initial parameters 
N = 128
A = 0.2
Δt = 0.01
tf = Float64(2000)
L = 2π;
k = 1;
h = 0;

using MAT
vars = matread(projectdir()*"/data/ClamondA2N128.mat")
X = vec(vars["X"])
Y = vec(vars["Y"]);
ϕ = vec(vars["F"])*sqrt(9.81);

# Looping parameters
algs = [Vern7(),Tsit5()];
smoothing = [true,false];
reltols = [1e-3,1e-4,1e-5,1e-6,1e-7,1e-8]


for alg ∈ algs
    for smoothes ∈ smoothing
        for reltol ∈ reltols

println(alg)
# Using Clamond Ak=0.2 Stokes wave 

if alg == Vern7()
    algname = "Vern7"
elseif alg == Tsit5()
    algname = "Tsit5"
end

description = "alg=$algname"*"_smooth=$smoothes"*"_reltol=$reltol"
# Create folders to save output
directory = projectdir()*"/plots/Longevolution/"*description
# Remove if already exists
rm(directory,recursive=true,force=true)
mkdir(directory)


# Run simulation
println("Starting run with $algname algorithm, with smoothing = $smoothes, and reltol = $reltol")
sol = @time runSim(N,X,Y,ϕ,Δt,tf,L,smoothing=smoothes,alg = alg,reltol=reltol)
data = Dict("solution "=> sol, "alg"=>alg,"smoothing"=>smoothes,"reltol"=>reltol)
@save directory*"/data.jld2" data

# Check energy + MWL
energy, MWL_check = computeEnergy(sol,N,Δt,tf)
gr()
energyplot = plot(0:Δt:tf,energy .- energy[1],xlabel=L"t \, (s)",ylabel =L"E - E_0 \quad (J/kg)",
legend=false,title= L"E_0: \quad "*@sprintf("%.4f ", energy[1])*L"(J/kg)",dpi=300)
savefig(energyplot,directory*"/energy.png")
MWLplot = plot(0:Δt:tf,MWL_check,xlabel=L"t \, (s)",ylabel ="MWL"*L" \quad (m)",legend=false,dpi=300)
savefig(MWLplot,directory*"/mwl.png")

# Spectra
using CubicSplines
using FFTW
gr()
# initial spectra
x = sol(0)[1:N];
y = sol(0)[N+1:2*N];
spline = CubicSpline(x,y);
ranges = range(minimum(x),maximum(x),length=2*N);
wavenum = fftfreq(2*N, 2*N /(2π))*2π;
surface = spline[ranges];
spectra1 = (abs.(fft(surface)/(2*N)).^2)[1:N] ;
spectra1[2:end] *= 2;

anim = @animate for i∈0:Δt*1000:tf
    x = sol(i)[1:N];
    y = sol(i)[N+1:2*N];
    spline = CubicSpline(x,y);
    ranges = range(minimum(x),maximum(x),length=2*N);
    wavenum = fftfreq(2*N, 2*N /(2π))*2π;
    surface = spline[ranges];
    spectra = (abs.(fft(surface)/(2*N)).^2)[1:N] ;
    spectra[2:end] *= 2;
    plot(wavenum[2:N],spectra[2:end],yscale=:log10,xscale=:log10,xaxis="cycles per 2π",yaxis = "η PSD (m^2 / cp2π)",label="current spectrum",
    title = @sprintf("Time: %.1f s", i),
    ylims=(1e-12,5e-2))
    plot!(wavenum[2:N],spectra1[2:end],label="initial spectrum")
end
gif(anim,directory*"/spectra.gif")
end
end
end

# Compare with Dold
doldData = load(projectdir()*"/plots/Longevolution/DoldA2N128MWL.jld2")
xD = doldData["x_o"];
yD = doldData["y_o"];
tD = doldData["t"];

x = xD[1,:];
y = yD[1,:];
spline = CubicSpline(x,y);
ranges = range(minimum(x),maximum(x),length=2*N);
wavenum = fftfreq(2*N, 2*N /(2π))*2π;
surface = spline[ranges];
spectra1 = (abs.(fft(surface)/(2*N)).^2)[1:N] ;
spectra1[2:end] *= 2;

function detrend(y)
    x = 1:length(y)
    coeffs = hcat(ones(length(x)), x) \ y
    return y .- (coeffs[1] .+ coeffs[2] .* x)
end
xdet = detrend(x)
ydet = detrend(y)
xSpectra = abs.(fft(xdet)/(N).^2)[2:N÷2]*2*N
plot(xSpectra,scale=:log10)
ySpectra = abs.(fft(ydet)/(N).^2)[2:N÷2]*2*N
plot(ySpectra[2:end],scale=:log10)


anim = @animate for i∈1:100:length(tD)
    x = xD[i,:];
    y = yD[i,:];
    spline = CubicSpline(x,y);
    ranges = range(minimum(x),maximum(x),length=2*N);
    wavenum = fftfreq(2*N, 2*N /(2π))*2π;
    surface = spline[ranges];
    spectra = (abs.(fft(surface)/(2*N)).^2)[1:N] ;
    spectra[2:end] *= 2;
    plot(wavenum[2:N],spectra[2:end],yscale=:log10,xscale=:log10,xaxis="cycles per 2π",yaxis = "η PSD (m^2 / cp2π)",label="current spectrum",
    title = @sprintf("Time: %.1f s", tD[i]),
    ylims=(1e-12,5e-2))
    plot!(wavenum[2:N],spectra1[2:end],label="initial spectrum")
end
gif(anim,projectdir()*"/plots/Longevolution/Doldspectra.gif")


anim = @animate for i∈500:600
    x = xD[i,:];
    y = yD[i,:];
    # spline = CubicSpline(x,y);
    # ranges = range(minimum(x),maximum(x),length=2*N);
    # wavenum = fftfreq(2*N, 2*N /(2π))*2π;
    # surface = spline[ranges];
    # spectra = (abs.(fft(surface)/(2*N)).^2)[1:N] ;
    # spectra[2:end] *= 2;
    # plot(wavenum[2:N],spectra[2:end],yscale=:log10,xscale=:log10,xaxis="cycles per 2π",yaxis = "η PSD (m^2 / cp2π)",label="current spectrum",
    # title = @sprintf("Time: %.1f s", tD[i]),
    # ylims=(1e-12,5e-2))
    # plot!(wavenum[2:N],spectra1[2:end],label="initial spectrum")
    scatter(x,y,aspect_ratio=1)
end
gif(anim,projectdir()*"/plots/Longevolution/Doldseries.gif")