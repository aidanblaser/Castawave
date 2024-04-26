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
reltols = [1e-9];
Nrange = [200];
Arange = 0.02:0.02:0.42

for N ∈ Nrange
    for A ∈ Arange 
        for smoothes ∈ smoothing
            for reltol ∈ reltols

                alg = Rodas4P(autodiff=false)
                if alg == Vern7()
                    algname = "Vern7"
                elseif alg == AB4()
                    algname = "AB4"
                elseif alg == ImplicitMidpoint(autodiff=false)
                    algname = "ImplicitMidpoint"
                elseif alg == Rodas4P(autodiff=false)
                    algname = "Rodas4P"
                elseif alg == Kvaerno5(autodiff=false)
                    algname = "Kvaerno5"
                elseif alg == SSPRK22()
                    algname = "SSPRK22"
                elseif alg == TRBDF2(autodiff=false)
                    algname = "TRBDF2"
                elseif alg == KenCarp5(autodiff=false)
                    algname = "KenCarp5"
                elseif alg == SSPRK22()
                    algname = "SSPRK22"
                end

                kH2 = A

                description = "N=$N"*"_kH2=$kH2"*"_alg=$algname"*"_smooth=$smoothes"*"_reltol=$reltol"
                # Create folders to save output
                directory = projectdir()*"/data/sims/Arange_N=200_reltol=1e-9_alg=$algname/"*description
                # Remove if already exists
                rm(directory,recursive=true,force=true)
                mkdir(directory)

                # Derive initial conditions 
                X,Y,ϕ,c = getIC(Inf,kH2,N÷2);
                L = 2π;
                # Start timer 
                timer = TimerOutput()
                Δt = 0.0001
                tf = 5*(2π/c)
                # Run simulation
                println("Starting run with N=$N, kH2 =$kH2,  $algname algorithm, with smoothing = $smoothes, and reltol = $reltol")
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
end

# For non-smoothed, make a table of
Nrange = range(20,stop=1000,step=20);
kH2range = 0.02:0.02:0.42
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
N = 200;

for kH2 in kH2range
    alg = Rodas4P(autodiff=false);
            if alg == Vern7()
                algname = "Vern7"
            elseif alg == Tsit5()
                algname = "Tsit5"
            elseif alg == SSPRK22()
                algname = "SSPRK22"
            elseif alg == TRBDF2(autodiff=false)
                algname = "TRBDF2"
            elseif alg == Rodas4P(autodiff=false)
                algname = "Rodas4P"
            end

    # Next, read in data, and make table of information
    description = "N=$N"*"_kH2=$kH2"*"_alg=$algname"*"_smooth=true"*"_reltol=1.0e-9"
    directory = projectdir()*"/data/sims/Arange_N=$N"*"_reltol=1e-9_alg=$algname/"*description
    data = load(directory*"/data.jld2")["data"]
    c = data["c"];
    if haskey(data,"Castawavesol") #check if code ran
        sol = data["Castawavesol"]
        if maximum(sol.t) ≈ 5*(2π)/c #check to make sure code didn't quit early
            energy = data["energy"]
            MWL = data["MWL"]
            c_est = data["c_est"]
            push!(energyRanges,(mean(energy[end-1000:end])-mean(energy[1:1000]))/(energy[1])*100)
            push!(MWLRanges,(mean(MWL[end-1000:end])-mean(MWL[1:1000])))
            push!(c_Ranges, (mean(c_est[end-1000:end])-c)/c * 100)
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
plottitle = "N = 200, tol=1e-9, alg = Rodas4P, smooth"
eplot = plot(kH2range[stableRange],energyRanges[stableRange],xlabel = "Wave slope",ylabel = "ΔE over five periods (%)",label="Castawave",title=plottitle)
plot!(kH2range[stableRange],energyDold[stableRange],label="Dold")
mwlplot = plot(kH2range[stableRange],MWLRanges[stableRange],xlabel="Wave slope",ylabel="ΔMWL over five periods (m)",label="Castawave",title=plottitle)
plot!(kH2range[stableRange],MWLDold[stableRange],label="Dold")
cplot = plot(kH2range[stableRange],c_Ranges[stableRange],xlabel="Wave slope",ylabel="Δc over five wave periods (%)",label="Castawave",title=plottitle)
plot!(kH2range[stableRange],cDold[stableRange],label="Dold")
khplot = plot(kH2range[stableRange],Waveheight_Ranges[stableRange],xlabel="Wave slope",ylabel="ΔH over five wave periods (%)",legend=false,title=plottitle)
plot!(kH2range[stableRange],kHDold[stableRange],label="Dold")
directory = projectdir()*"/data/sims/Arange_N=200_reltol=1e-9_alg=Rodas4P/"
savefig(eplot, directory*"eplot.png")
savefig(mwlplot, directory*"mwlplot.png")
savefig(cplot, directory*"cplot.png")
savefig(khplot, directory*"khplot.png")

# Save for reltol 1e-7
stableRange7 = stableRange
energyRange7 = energyRanges
MWLRanges7 = MWLRanges
c_Ranges7 = c_Ranges

# Save for reltol 1e-5
stableRange5 = stableRange
energyRange5 = energyRanges 
MWLRanges5 = MWLRanges
c_Ranges5 = c_Ranges

# Save for reltol 1e-6
plot(Nrange,energyRanges,label="Castawave 1e-6 tol",ylabel="ΔEnergy over 5 periods")
plot!(Nrange,energyRange7,label="Castawave 1e-7 tol")
plot!(Nrange[end-2],coalesce.(energyRange5, missing),label="Castawave 1e-5 tol")


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
