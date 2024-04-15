using DrWatson
@quickactivate "Castawave"
using Test
using TimerOutputs

# Import solver and Clamond IC function 
include(projectdir()*"/src/MainSolver.jl")
include(projectdir()*"/src/ClamondIC.jl")
include(projectdir()*"/src/DoldMono.jl")


# First test: run through various numbers of points for 5 wave periods at kH/2 = 0.1
Nrange = range(20,stop=1000,step=20);
smoothing=[true,false];
reltols = [1e-6];

for N ∈ Nrange
    for smoothes ∈ smoothing
        for reltol ∈ reltols

            alg = Vern7();
            if alg == Vern7()
                algname = "Vern7"
            elseif alg == Tsit5()
                algname = "Tsit5"
            end

            kH2 = 0.1

            description = "N=$N"*"_alg=$algname"*"_kH2=$kH2"*"_smooth=$smoothes"*"_reltol=$reltol"
            # Create folders to save output
            directory = projectdir()*"/data/sims/VaryPoints-Ak=0.1/"*description
            # Remove if already exists
            rm(directory,recursive=true,force=true)
            mkdir(directory)

            # Derive initial conditions 
            X,Y,ϕ,c = getIC(Inf,kH2,N÷2);
            L = 2π;
            # Start timer 
            timer = TimerOutput()
            Δt = 0.01
            tf = 5*(2π/c)
            # Run simulation
            println("Starting run with $algname algorithm, with smoothing = $smoothes, and reltol = $reltol")
            # If error, go to next datapt
            try
            sol = @timeit timer "simulation" runSim(N,X,Y,ϕ,Δt,tf,L,smoothing=smoothes,alg = alg,reltol=reltol)
            # Check energy + MWL
            energy, MWL_check, phasespd = computeEnergy(sol,N,Δt,tf,L)
            gr()
            energyplot = plot(0:Δt:tf,energy .- energy[1],xlabel=L"t \, (s)",ylabel =L"E - E_0 \quad (J/kg)",
            legend=false,title= L"E_0: \quad "*@sprintf("%.4f ", energy[1])*L"(J/kg)",dpi=300)
            savefig(energyplot,directory*"/energy.png")
            MWLplot = plot(0:Δt:tf,MWL_check,xlabel=L"t \, (s)",ylabel ="MWL"*L" \quad (m)",legend=false,dpi=300)
            savefig(MWLplot,directory*"/mwl.png")
            cplot = plot(0:Δt:tf,(phasespd .- c)./c .*100,xlabel=L"t \, (s)",ylabel = "Δc (%)",legend=false,dpi=300)
            savefig(cplot,directory*"/phasespd.png")

            # Run same simulation in Dold
            # Have to modify ww_write.f to have correct N 
            # Read the Fortran file
            file_content = read(projectdir()*"/src/dold_files/ww_write.f",String)
            # Search for (niv = *number*)
            pattern = r"niv\s*=\s*\d+"
            match_content = match(pattern, file_content)
            content_inside_parentheses = match_content.match
            # Replace it with new value for N 
            modified_string = replace(file_content, content_inside_parentheses => "niv = $N")
            write(projectdir()*"/src/dold_files/ww_write.f", modified_string)
            xD,yD,ϕD,t = monoSim(X,Y,ϕ,tf);

            # Save data
            data = Dict("Castawavesol"=> sol,"Timer"=>timer,"N"=>N,"kH2"=>kH2, "alg"=>alg,
            "smoothing"=>smoothes,"reltol"=>reltol,"energy"=>energy,"MWL"=>MWL_check,"c_est"=>phasespd,"c"=>c,
            "xDold"=>xD,"yDold"=>yD,"ϕDold"=>ϕD,"tDold"=>t)
            @save directory*"/data.jld2" data
            catch
            println("Solution unstable")
            # Run same simulation in Dold
            # Have to modify ww_write.f to have correct N 
            # Read the Fortran file
            file_content = read(projectdir()*"/src/dold_files/ww_write.f",String)
            # Search for (niv = *number*)
            pattern = r"niv\s*=\s*\d+"
            match_content = match(pattern, file_content)
            content_inside_parentheses = match_content.match
            # Replace it with new value for N 
            modified_string = replace(file_content, content_inside_parentheses => "niv = $N")
            write(projectdir()*"/src/dold_files/ww_write.f", modified_string)
            xD,yD,ϕD,t = monoSim(X,Y,ϕ,tf);
            data = Dict("N"=>N,"kH2"=>kH2, "alg"=>alg,
                "smoothing"=>smoothes,"reltol"=>reltol,"c"=>c,
                "xDold"=>xD,"yDold"=>yD,"ϕDold"=>ϕD,"tDold"=>t)
            @save directory*"/data.jld2" data
            end

        end
    end
end

# For non-smoothed, make a table of
Nrange = range(20,stop=1000,step=20);
energyRanges = [];
MWLRanges = [];
c_Ranges = [];
runtime_Ranges = [];
unstable_Ranges = [];
Waveheight_Ranges = [];
MWLDold = [];
energyDold = [];
cDold = [];
kHDold = [];

for N in Nrange
    alg = Vern7();
            if alg == Vern7()
                algname = "Vern7"
            elseif alg == Tsit5()
                algname = "Tsit5"
            end

    # Next, read in data, and make table of information
    description = "N=$N"*"_alg=$algname"*"_kH2=0.1"*"_smooth=false"*"_reltol=1.0e-6"
    directory = projectdir()*"/data/sims/VaryPoints-Ak=0.1/"*description
    data = load(directory*"/data.jld2")["data"]
    c = data["c"];
    if haskey(data,"Castawavesol") #check if code ran
        sol = data["Castawavesol"]
        if maximum(sol.t) ≈ 5*(2π)/c #check to make sure code didn't quit early
            energy = data["energy"]
            MWL = data["MWL"]
            c_est = data["c_est"]
            push!(energyRanges,(energy[end]-energy[1])/(energy[1])*100)
            push!(MWLRanges,(MWL[end]-MWL[1]))
            push!(c_Ranges, (c_est[end]-c)/c * 100)
            push!(unstable_Ranges,0.0)
            push!(Waveheight_Ranges,((maximum(sol(5*(2π)/c)[N+1:2*N]) - minimum(sol(5*(2π)/c)[N+1:2*N]))- (maximum(sol(0)[N+1:2*N]) - minimum(sol(0)[N+1:2*N])))/((maximum(sol(0)[N+1:2*N]) - minimum(sol(0)[N+1:2*N])))*100)
        else
            push!(energyRanges,NaN)
            push!(MWLRanges,NaN)
            push!(c_Ranges,NaN)
            push!(unstable_Ranges,maximum(sol.t))
            push!(Waveheight_Ranges,NaN)
        end
    else
        push!(energyRanges,NaN)
        push!(MWLRanges,NaN)
        push!(c_Ranges,NaN)
        push!(unstable_Ranges,NaN)
        push!(Waveheight_Ranges,NaN)
    end
    # Load dold solution 
    xD = data["xDold"]
    yD = data["yDold"]
    ϕD = data["ϕDold"]
    tD = data["tDold"]
    # Get MWL, energy, deviations
    xDend = xD[end,:];
    xDstart = xD[1,:];
    xDξstart = DDI1(xDstart,N,2π,1)
    xDξend = DDI1(xDend,N,2π,1)
    yDξstart = DDI1(yD[1,:],N,0,1)
    yDξend = DDI1(yD[end,:],N,0,1)
    MWLend = sum(yD[end,:].*xDξend)./N
    MWLstart = sum(yD[1,:].*xDξstart)./N
    push!(MWLDold,MWLend-MWLstart)
    kH = (maximum(yD[end,:]) - minimum(yD[end,:]) - maximum(yD[1,:]) + minimum(yD[1,:]))/(maximum(yD[1,:]) - minimum(yD[1,:])) * 100;
    push!(kHDold, kH)
    # Get velocities
    dXstart, dYstart, dϕstart = fixedTimeOperations(N,xD[1,:],yD[1,:],ϕD[1,:],2π,0.0,false)
    dXend, dYend, dϕend = fixedTimeOperations(N,xD[end,:],yD[end,:],ϕD[end,:],2π,0.0,false)
    Bξstart = dXstart .* yDξstart .- dYstart .* xDξstart
    Bξend = dXend .* yDξend .- dYend .* xDξend
    integrand_start = -1/2 * ϕD[1,:] .* Bξstart
    integrand_end = -1/2 *ϕD[end,:].*Bξend
    KEstart = sum(integrand_start)*(2π)/N
    KEend = sum(integrand_end)*(2π)/N
    PEstart = sum(GRAVITY/2 * (yD[1,:]).^2 .* xDξstart)*2π/N
    PEend = sum(GRAVITY/2 * (yD[end,:]).^2 .* xDξend)*2π/N
    push!(energyDold, (KEend + PEend - KEstart - PEstart)/(KEstart - PEstart)*100)
    push!(cDold, median((Bξend ./ yDξend .- Bξstart ./ yDξstart)./(Bξstart ./ yDξstart)*100))


end
# Find ranges where solution works

gr()
stableRange = findall(isequal(0.0),unstable_Ranges)
plottitle = "kH/2 = 0.1, reltol=1e-6, smoothed = false"
eplot = plot(Nrange[stableRange],energyRanges[stableRange],xlabel = "N",ylabel = "ΔE over five periods (%)",label="Castawave",title=plottitle)
plot!(Nrange[stableRange],energyDold[stableRange],label="Dold")
mwlplot = plot(Nrange[stableRange],MWLRanges[stableRange],xlabel="N",ylabel="ΔMWL over five periods (m)",label="Castawave",title=plottitle)
plot!(Nrange[stableRange],MWLDold[stableRange],label="Dold")
cplot = plot(Nrange[stableRange],c_Ranges[stableRange],xlabel="N",ylabel="Δc over five wave periods (%)",label="Castawave",title=plottitle)
plot!(Nrange[stableRange],cDold[stableRange],label="Dold")
khplot = plot(Nrange[stableRange],Waveheight_Ranges[stableRange],xlabel="N",ylabel="ΔH over five wave periods (%)",legend=false,title=plottitle)
plot!(Nrange[stableRange],kHDold[stableRange],label="Dold")
savefig(eplot, projectdir()*"/data/sims/VaryPoints-Ak=0.1/eplot.png")
savefig(mwlplot, projectdir()*"/data/sims/VaryPoints-Ak=0.1/mwlplot.png")
savefig(cplot, projectdir()*"/data/sims/VaryPoints-Ak=0.1/cplot.png")
savefig(khplot, projectdir()*"/data/sims/VaryPoints-Ak=0.1/khplot.png")


description = "N=160"*"_alg=Vern7"*"_kH2=0.3"*"_smooth=false"*"_reltol=1.0e-6"
    directory = projectdir()*"/data/sims/VaryPoints/"*description
    data = load(directory*"/data.jld2")["data"]
sol = data["Castawavesol"]
scatter(sol(9.5)[1:160],sol(9.5)[161:320])
MWL = data["MWL"]
plot(MWL)
MWL[end]

t = 0.5:0.001:0.57169
using LaTeXStrings
using Printf
gr()
anim = @animate for i ∈ t
    xvals = mod.(sol(i)[1:N],2π)
    yvals = sol(i)[N+1:2*N]
    ϕvals = sol(i)[2*N+1:3*N]

    scatter(xvals,yvals, legend = false,
        framestyle= :box,background_color="black", markerstrokewidth=0, markersize=1,
        dpi = 300, xlabel=L"x \,(m)",ylabel=L"z \,(m)", title= @sprintf("Time: %.4f s", i),
        xlims=(1.5,2),ylims = (0,0.7))
    #plot!([maxVec[Int(t÷Δt + 1)],maxVec[Int(t÷Δt + 1)]],[-2,2],linewidth=3)
end every 1
gif(anim, projectdir()*"/plots/animation.gif", fps=10)
