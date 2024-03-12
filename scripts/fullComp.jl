using DrWatson
@quickactivate "Castawave"
using LaTeXStrings
using Printf
using JLD2

# reads in all files, and compares
castawave = load(projectdir()*"/data/CastawaveN512.jld2")
dold = load(projectdir()*"/data/DoldN512.jld2")

# Start with 1st simulation for now
simNumber = 6

doldX = dold["XArray"][simNumber];
doldY = dold["YArray"][simNumber];
doldϕ = dold["ϕArray"][simNumber];
doldt = dold["tArray"][simNumber];

T,N = size(doldX);

castawave1st = castawave["SolArray"][simNumber];

# create arrays from castawave on same grid
rmseXArray = [];
rmseYArray = [];
rmseϕArray = [];
for i in 1:T
    castawaveX = castawave1st(doldt[i])[1:N];
    castawaveY = castawave1st(doldt[i])[N+1:2*N];
    castawaveϕ = castawave1st(doldt[i])[2*N+1:end];

    rmseX = sqrt(sum((doldX[i,:] .- castawaveX).^2)/N);
    rmseY = sqrt(sum((doldY[i,:] .- castawaveY).^2)/N);
    rmseϕ = sqrt(sum((doldϕ[i,:] .- castawaveϕ).^2)/N);

    append!(rmseXArray,rmseX)
    append!(rmseYArray,rmseY)
    append!(rmseϕArray,rmseϕ)
end

plot(doldt,rmseXArray,xlabel="time",ylabel="rmseX (m)")
plot(doldt,rmseYArray,xlabel="time",ylabel="rmseY (m)")
plot(doldt,rmseϕArray)

simpsons_rule_periodic(doldX[end,:],doldY[end,:])
simpsons_rule_periodic(doldX[1,:],doldY[1,:])

doldspectra = load(projectdir()*"/data/doldspectra.jld2")
spectraDold = doldspectra["spectra"]
time = doldspectra["t"]

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
Y = [(cos(k * α*L / n )) * A + 1/6*A^4*k^3*cos(k*α*L/n) .+ 0*(A^2*k / 2).+ 0*A^4*k^3 * 1/2 for α in 1:n]
ϕ = [sqrt(GRAVITY/k) * A * exp.(k*Y[α]) * sin(k*X[α]) for α in 1:n]

sol= runSim(n, X, Y, ϕ, Δt, Float64(tf),L,h)
using CubicSplines
using FFTW
spectra = [];
kranges = [];
for t∈ time
    yvals = sol(t)[n+1:2*n]
    xvals = sol(t)[1:n]
    interp = CubicSpline(xvals,yvals)
    xrange = range(xvals[1],stop = xvals[end],length=1000)
    eta = interp[xrange]
    FFT = fft(eta).*(xrange[2]-xrange[1])
    krange = fftfreq(length(xrange),xrange[2]-xrange[1])
    push!(spectra,abs.(FFT[2:end÷2]).^2)
    push!(kranges,krange[2:end÷2])
end
anim = @animate for i ∈ 1:length(spectra)
    times = @sprintf("%0.1f",time[i]);
    plot(spectra[i],xaxis=:log,yaxis=:log,xlabel = "cycles per 2π meters",ylabel=L"|\hat{\eta}|^2 \quad (m)",
    title = "Time = $times (s)",ylims = (1e-14,1e0),label="Castawave Spectra A=0.3")
    plot!(spectraDold[i],label="Dold Spectra A=0.3")
end 
gif(anim, "Plots/spectra.gif", fps = 10)
