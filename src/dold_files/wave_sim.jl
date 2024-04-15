using DrWatson
@quickactivate "DrWatson Example"
using FFTW
using Plots
using DelimitedFiles

#= Wave breaking criteria (Aidan's Tests)
#--------------------------------------
# June 2021. Nick Pizzo
# Overview: Our goal is to examine geometric, kinematic and dynamic
# properties of focusing and breaking waves. To do this we 
# employ Dold's code which is a mixed Eulerian-Lagrangian solver
# which can model waves up to the point of surface reconnection. See 
# Dold (1992)
#---------------------------------------
# We wish to resolve: 
# Geometry: Wave height, slope, curvature 
# Kinematics: particle velocity, acceleration, phase speed 
# Dynamics: Kinetic, potential energy
# Kinematic breaking criteria (ratio of particle speed to phase speed)
# approaching breaking at highest & steepest part of wave
#--------------------------------------- 
# Simulation parameters
=#

N = 512; # Number of surface points
ML = 100; # Physical length of channel
xₒ = -30; # Phase shift to match initial conditions with lab data
x =collect( xₒ : ML/N : xₒ + (ML) * (1 - 1/N)); # domain
wl = ML; # wl parameter is used in Dold code
g = 9.81; # acc due to gravity. make sure Dold is set to same value
k = (2π)^2 / g; # central wavenumber
w = sqrt(g*k); # associated central angular freq - deep water
cg = g / 2 / w; # group velocity 
D = 1; # bandwidth
xb = 20/D; #focusing location
tb = xb/cg; # focusing time
S = 0.05; # linear prediction max slope at focusing
M = 32*4; # number of modes
I = im; # define imaginary number
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
eta = zeros(1, length(x)); 
for j ∈ 1 : length(x)
    eta[1,j] = sum( real( 2 .* out[2:end].*
    exp.( I .* ( 2 * π .* freq[2:end].*(x[j] .+ 30) ) )));
end
# check S for the windowed wave form
eps2 = sum((abs.(2 .*out[2:end]) .* 2 .* π .* freq[2:end] )); #
# now find phi
phi = zeros(1, length(x)); 
for j ∈ 1 : length(x)
    phi[1,j] = sum(imag.(2 .* sqrt(g) .* out[2:end] ./ 
        sqrt.(2 .* π .* freq[2:end]).*
        exp.( I .* (( 2 .* π .* freq[2:end] ).*
        (x[j] - xₒ) ))));
end
# BW=powerbw(eta,2*pi*Fs,[],10)/k ; # check bandwidth

x_f = x .- xₒ; # shift x-axis so it goes from 0 to wl
y_f = eta;
f_f = phi;

plot( x_f, [y_f',f_f'], xlabel='x', label=["η" "ϕ"],
        legendfontsize=18)
#set(gca, 'fontsize', 22)
#l1 = legend('$\eta$' , '$\phi$');
#set(l1, 'interpreter' , 'latex')


## 
BW = 1; # this is set to 1 for now
tl = (xb+20)*D;
#tl = (xb+20)*D; # length of simulation

# Make sure you're in right directory
cd(projectdir()) 
cd("scripts")

rm("xc.txt")
rm("yc.txt")
rm("fc.txt")
rm("bw.txt")
rm("wl.txt")
rm("S.txt")
rm("bw0.txt")
rm("S0.txt")
# delete k.txt
writedlm("xc.txt",x_f)
writedlm("yc.txt", y_f)
writedlm("fc.txt", f_f)
writedlm("S.txt", eps2)
writedlm("S0.txt", S)
writedlm("bw.txt", BW)
writedlm("bw0.txt", BW)
writedlm("C.txt", BW)
writedlm("wl.txt", wl)
writedlm("tl.txt", tl)
# save k.txt k -ascii
run(`./run2.sh`);

# Open files and save as variables
x_o = readdlm("x.txt");
y_o = readdlm("y.txt");
p_o = readdlm("phi.txt");
to = readdlm("t.txt");
t = to[:,1];


# re-dimensionalize variables
L = wl / 2 / π; #length scale
T = sqrt( wl / 2 / π ); # time scale
xi = 0 : wl / N : wl - wl / N; # spatial grid
# look at the data

@gif for i ∈ 1:2:length(t)
    time = round(T.*t[i],digits=3);
    plot(x_o[i,:].*L , y_o[i,:].*L, xlims=(0,wl),
    ylims = (-0.015,0.015), xlabel="x (m)",ylabel="η (m)",
    title = "S = $S , Δ = $D, t = $time s")
end 
    #=for i ∈ 1 : 5 : length(t)
plot( x_o(:,i) * L , y_o(:,i) * L , 'k')
axis([0 wl -0.01 0.01])
grid on 
title(sprintf('time = %0.2f', T * t(i) ));
xlabel('X (m)') 
ylabel('Y (m)')
set(gca , 'fontsize' , 28)
hold off 
pause(0.1)
end
=#