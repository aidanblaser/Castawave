using DrWatson
@quickactivate "Castawave"
using LaTeXStrings
using Printf
using JLD2

# reads in all files, and compares
castawave = load(projectdir()*"/data/CastawaveN128.jld2")
dold = load(projectdir()*"/data/DoldN128.jld2")

# Start with 1st simulation for now
simNumber = 8

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
