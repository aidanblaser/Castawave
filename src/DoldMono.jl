# General Dold Simulation for Monochromatic
# Differs from DoldSim as it mostly uses single waves

using DrWatson
@quickactivate "Castawave"
using FFTW
using Plots
using DelimitedFiles
using Interpolations
using LaTeXStrings
using Printf
using Statistics

function monoSim(N,X,Y,ϕ,T,plotting = nothing)
directory = pwd();
# Adjustable Parameters
g = 9.81; # acc due to gravity. make sure Dold is set to same value
k = 1;
tl = T ; # length of simulation (seconds)
ML = 2π*k; # Physical length of channel
x =collect( 0 : ML/N : (ML) * (1 - 1/N)); # domain
BW = 0; # bandwidth (set to 0 for no focusing packet)
wl = ML; # wl parameter is used in Dold code

x_f =X;
y_f = Y;
f_f = ϕ;

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

cd(projectdir()) 
#look at the data
if plotting !=nothing
anim = @animate for i ∈ 1:length(t)
    time = @sprintf("%0.1f",t[i]);
    plot(x_o[i,:] , y_o[i,:],  xlabel="x (m)",ylabel="η (m)",
    title = "S = $S, t = $time s",label="Numerics",
    xlims=(minimum(x_o),maximum(x_o)),
    ylims=(minimum(y_o.-1.5.*v),maximum(y_o.+1.5*v)),
    framestyle= :box,background_color="black",legend=false,)
    #plot!(x_o[1:i,250].*L,y_o[1:i,250].*L,label="",color="red")
    #scatter!([x_o[i,250].*L],[y_o[i,250].*L],label="")
end 
gif(anim, "Plots/monoSim.gif", fps = 6)
end

cd(directory)

return x,y,ϕ,t;
end

# N = 128; #number of points
# A = 0.2;
# k = 1;
# g = 9.81;
# # X = collect(range(0,step=2π/N,length=N));
# # Y = A*cos.(k*X) .+ A^3 * k^2 / 8 * cos.(k*X) .+ A^2*k/2 * cos.(2*k*X) .+ 3 *A^3 *k^2 / 8 * cos.(3*k*X);
# # ϕ = A*sqrt(g/k)*sin.(k*X);

# X = [2π*x/ N - A*sin.(2π*x/N) for x∈ 1:N ];
# Y = [A*cos.(2π*x/N) for x∈ 1:N];
# ϕ = [sqrt(g)*A*exp(Y[i])*sin.(X[i]) for i∈ 1:N ];

# using MAT
# data = matread(projectdir()*"/Data/ClamondA2N128.mat");
# X = vec(data["X"]);
# Y = vec(data["Y"]);
# ϕ = vec(data["F"])*sqrt(9.81);



# x_o,y_o,t = @time monoSim(X,Y,ϕ,100.)
# using JLD2
# jldsave(projectdir()*"/Data/DoldA2N128MWL.jld2";x_o,y_o,t)



# cd(projectdir)
# # Plot results
# anim = @animate for i ∈ 1:length(t)
#     time = @sprintf("%0.1f",t[i]);
#     plot(x_o[i,:] , y_o[i,:],  xlabel="x (m)",ylabel="η (m)",
#     title = "A = $A, t = $time s",label="Numerics",
#     xlims=(minimum(x_o),maximum(x_o)),
#     ylims=(minimum(y_o),maximum(y_o)),
#     framestyle= :box,background_color="black",legend=false,)
# end 
# gif(anim, "Plots/monoSim.gif", fps = 6)

# # Interpolate to get η(x)
# using CubicSplines
# using FFTW
# spectra = [];
# kranges = [];
# for i ∈ 1:length(t)
#     interp = CubicSpline(x_o[i,:],y_o[i,:])
#     xrange = range(x_o[i,1],stop = x_o[i,end],length=1000)
#     eta = interp[xrange]
#     FFT = fft(eta).*(xrange[2]-xrange[1])
#     krange = fftfreq(length(xrange),xrange[2]-xrange[1])
#     push!(spectra,abs.(FFT[2:end÷2]).^2)
#     push!(kranges,krange[2:end÷2])
# end
# anim = @animate for i ∈ 1:length(t)
#     time = @sprintf("%0.1f",t[i]);
#     plot(spectra[i],xaxis=:log,yaxis=:log,xlabel = "cycles per 2π meters",ylabel=L"|\hat{\eta}|^2 \quad (m)",
#     title = "Time = $time (s)",ylims = (1e-14,1e0),label="Dold Spectra A=0.3")
# end 
# gif(anim, "Plots/spectra.gif", fps = 30)

# jldsave(projectdir()*"/data/DoldSpectra.jld2"; spectra,t)

# X = [2π*x/ N - 0.01*sin.(2π*x/N) for x∈ 1:N ];
# Y = [0.01*cos.(2π*x/N) for x∈ 1:N];
# ϕ = [sqrt(g)*0.01*sin.(2π*x/N) for x∈ 1:N ];

# plot(t,mean(y_o,dims=2))
# using Trapz
# mwl = zeros(length(t));
# for i ∈ 1:length(t)
#     mwl[i] = trapz(x_o[i,:],y_o[i,:]/(2π));
# end

# plot(t,mwl)

# plot(diff(t)[end-20:end])
