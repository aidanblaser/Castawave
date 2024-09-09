using DrWatson
@quickactivate "Castawave"
using Plots
using FFTW

N = 2048; # Number of surface points (must be changed in ww_write.f)
ML = 100; # Physical length of channel
xₒ = -30; # Phase shift to match initial conditions with lab data
x =collect( xₒ : ML/N : xₒ + (ML) * (1 - 1/N)); # domain
wl = ML; # wl parameter is used in Dold code
g = 9.81; # acc due to gravity. make sure Dold is set to same value
k = (2π)^2 / g; # central wavenumber
w = sqrt(g*k); # associated central angular freq - deep water
cg = g / 2 / w; # group velocity 
D = 0.8; # bandwidth
xb = 20/D; #focusing location
tb = xb/cg; # focusing time
S = 0.35; # linear prediction max slope at focusing
M = 32*4; # number of modes
kn = zeros(M,1); # allocate space for wavenumbers
wn = zeros(M,1); # allocate space for frequencies
an = zeros(M,1); # allocate space for amplitudes
for i ∈ 1 : M
    kn[i] = k * ( 1 + D * (i - M/2) / M ); # define wavenumbers
    wn[i] = sqrt( g * kn[i] ); # define associated freq
    an[i] = S / kn[i] / M; # define amplitudes
end
# Define surface displacement eta. Linear dispersive focusing 
eta1 = zeros(length(x),1);
for i ∈ 1 : M
    eta1[:,1] += an[i] * cos.( kn[i] .* (x .- xb) .+ wn[i] .* tb );
end

# now window the initial waveform so we have 1 wave group
# This will be updated, it's proportional to delta k. 
filt1 = 1/2 * (tanh.(0.25*(x.- -15))-tanh.(0.25.*(x.-10))); 
eta2 = eta1.*filt1;
# 
# ##
# hold on
# plot(x, eta2)
# xlabel('x','interpreter','latex')
# ylabel('$\eta$','interpreter','latex')
# l1 = legend('Initial wave form','Windowed wave form');
# set(l1,'interpreter', 'latex')
# set(gca,'fontsize',22)
# xlim([min(x) max(x)])
# ##

## we now want to find the velocity potential
Fs = 1 / ( abs(x[2] - x[1]) ); # define sampling frequency in space
function positiveFFT(x,Fs)
    N=length(x); 
    k=collect(0:N-1); 
    T=N/Fs; 
    freq=k/T; #create the frequency range 
    X=FFTW.fft(x)/N; # normalize the data
    cutOff = Int(ceil(N/2)); 
    X = X[1:cutOff]; 
    freq = freq[1:cutOff];
    return X, freq
end
out,freq = positiveFFT(eta2, Fs); #perform fft
# check fft on eta
eta = zeros(length(x)); 
for j ∈ 1 : length(x)
    eta[j] = sum( real( 2 .* out[2:end].*
    exp.( im .* ( 2 * π .* freq[2:end].*(x[j] .+ 30) ) )));
end
# check S for the windowed wave form
eps2 = sum((abs.(2 .*out[2:end]) .* 2 .* π .* freq[2:end] )); #
# now find phi
phi = zeros(length(x)); 
for j ∈ 1 : length(x)
    phi[j] = sum(imag.(2 .* sqrt(g) .* out[2:end] ./ 
        sqrt.(2 .* π .* freq[2:end]).*
        exp.( im .* (( 2 .* π .* freq[2:end] ).*
        (x[j] - xₒ) ))));
end
# BW=powerbw(eta,2*pi*Fs,[],10)/k ; # check bandwidth

x_f = x .- xₒ; # shift x-axis so it goes from 0 to wl
y_f = eta;
f_f = phi;

plot(x_f,y_f,label="Positions")
plot!(x_f,f_f,label="ϕ")


include(projectdir()*"/src/MainSolver.jl")
include(projectdir()*"/src/ClamondIC.jl")
include(projectdir()*"/src/DoldMono.jl")

tf = 4.28
L = 100.0

L̃ = L / 2π

X = x_f ./ L̃;
Y = y_f ./ L̃;
ϕ = 9.81^(3/2) .* f_f ./ (L̃^3/2);

h = 0.0
Δt = 1e-3
Xfull, Yfull, ϕfull, t = @time runSim(N, X, Y, ϕ, Δt, Float64(tf),2π,h,ϵ = tol,smoothing=smoothed)

gr()
function visualize(interval::Int, fps::Int)
    anim = @animate for i ∈ enumerate(t)
        plot(mod.(Xfull[i[1],:],L),Yfull[i[1],:], legend = false,
         framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1,
          dpi = 300, xlabel=L"x \,(m)",ylabel=L"z \,(m)", title= @sprintf("Time: %.1f s", i[2]),
          xlims=(2.5,3.1),ylims=(-0.007,0.007))
    end every interval
    gif(anim, projectdir()*"/plots/RK4noSmooth.gif", fps=fps)
end
visualize(10, 10)

include(projectdir()*"/src/DoldMono.jl")
xd, yd, ϕd, t =   monoSim(N,x_f,y_f,f_f,4.5)

scatter(xd[1,:],yd[1,:])
scatter(x_f,y_f)

directory = pwd();
# Adjustable Parameters
g = 9.81; # acc due to gravity. make sure Dold is set to same value
k = 1;
tl = 100. ; # length of simulation (seconds)
ML =100.; # Physical length of channel
x =collect( 0 : ML/N : (ML) * (1 - 1/N)); # domain
BW = 0; # bandwidth (set to 0 for no focusing packet)
wl = ML; # wl parameter is used in Dold code


# Make sure you're in right directory
cd(projectdir()) 
cd("src/dold_files")

rm("xc.txt",force=true)
rm("yc.txt",force=true)
rm("fc.txt",force=true)
rm("bw.txt",force=true)
rm("wl.txt",force=true)
rm("S.txt",force=true)
rm("bw0.txt",force=true)
rm("S0.txt",force=true)
rm("C.txt",force=true)
rm("tl.txt",force=true)
# delete k.txt
writedlm("xc.txt",x_f)
writedlm("yc.txt", y_f)
writedlm("fc.txt", f_f)
#writedlm("bw.txt", BW)
#writedlm("bw0.txt", BW)
#writedlm("C.txt", BW)
writedlm("wl.txt", wl)
writedlm("tl.txt", tl)
writedlm("S.txt", 0.05)
writedlm("S0.txt", 0.05)
# save k.txt k -ascii
run(`./run2.sh`);

# Open files and save as variables
x = readdlm("x.txt");
y = readdlm("y.txt");
ϕ = readdlm("phi.txt");
to = readdlm("t.txt");
t = to[:,1];
plotlyjs()
plot(x[end,:],y[end,:])
t[end]
visualize(10, 10)


jldsave(projectdir()*"/data/packetbreakN2048.jld2",x = x,y= y,ϕ= ϕ,t= t,N= N,S=S) 
using HDF5
file = h5open(projectdir()*"/data/packetbreakN2048.h5","w")
file["x"] = x 
file["y"] = y 
file["ϕ"] = ϕ
file["t"] = t 
file["N"] = N 
close(file)