using DrWatson
@quickactivate "Castawave"
using FFTW
using DelimitedFiles
using LinearAlgebra
using DSP
using Printf


include(projectdir()*"/src/HelperFunctions.jl")

function DoldSim(X,Y,ϕ,T,L)
    #=
    DoldSim is the main simulation function which inputs initial X,Y, and ϕ values
    and simulates its evolution until time T

    Inputs:
    X (Vector): Vector of initial X values of each Lagrangian particle (meters)
    Y (Vector): Vector of initial Y values of each Lagrangian particle (meters)
    ϕ (Vector): Vector of initial ϕ values of each Lagrangian particle (velocity potential)
    T (Float64): Time of simulation in seconds
    L (Float64): Length of periodic domain in meters

    Outputs:
    x (Matrix): Output matrix of X vector at each time  
    y (Matrix): Output matrix of Y vector at each time 
    ϕ (Matrix): Output matrix of ϕ vector at each time
    t (Vector): Output vector of timesteps where data is printed
    =#
    directory = pwd();
    # Adjustable Parameters
    g = 9.81; # acc due to gravity. make sure Dold is set to same value
    N = length(X); # Number of surface points
    wl = L # length of domain 
    println(projectdir())
    # Change ww_write.f if needed for correct N               
    # Read the Fortran file
    file_content = read(projectdir()*"/src/dold_files/ww_write_overturn.f",String)
    # Search for (niv = *number*)
    pattern = r"niv\s*=\s*\d+"
    match_content = match(pattern, file_content)
    content_inside_parentheses = match_content.match
    # Replace it with new value for N 
    modified_string = replace(file_content, content_inside_parentheses => "niv = $N")
    write(projectdir()*"/src/dold_files/ww_write_overturn.f", modified_string)


    x_f = X;
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
    rm("u.txt",force=true)
    rm("v.txt",force=true)
    rm("x.txt",force=true)
    rm("y.txt",force=true)
    rm("phi.txt",force=true)
    rm("enp.txt",force=true)
    rm("enk.txt",force=true)
    # delete k.txt
    writedlm("xc.txt",x_f)
    writedlm("yc.txt", y_f)
    writedlm("fc.txt", f_f)
    #writedlm("bw.txt", BW)
    #writedlm("bw0.txt", BW)
    #writedlm("C.txt", BW)
    writedlm("wl.txt", wl)
    writedlm("tl.txt", T)
    writedlm("S.txt", 0.1)
    writedlm("S0.txt", 0.1)
    # save k.txt k -ascii
    run(`./run2overturn.sh`);

    # Open files and save as variables
    x = readdlm("x.txt");
    y = readdlm("y.txt");
    ϕ = readdlm("phi.txt");
    to = readdlm("t.txt");
    u = readdlm("u.txt");
    v = readdlm("v.txt");
    KE = vec(readdlm("enk.txt"));
    PE = vec(readdlm("enp.txt"));
    t = to[:,1];

    # rescale variables to physical units 
    # re-dimensionalize variables
    L̃ = wl / 2 / π; #length scale
    T̃ = sqrt( wl / 2 / π ); # time scale

    x *= L̃
    y *= L̃
    ϕ *= L̃^(3/2)
    t *= T̃
    u *= L̃/T̃;
    v *= L̃/T̃;
    KE *= (L̃/T̃)^2;
    PE *= (L̃/T̃)^2;

    #return to orignal directory
    cd(directory)

    return x,y,ϕ,t,u,v,KE,PE;
end

function getIC(kd,kH2,N;tol=1e-14)
    #=
    Wrapper function for SSGW, which works in dimensionless units.
    getIC computes initial X, Y, ϕ values to high precision for an irrotational Stokes wave

    Inputs:
    kd (Float64): Depth of fluid (set to inf if deep-water)
    kH2 (Float64): Half the distance between wave crest and trough normalized by k (similar to slope)
    N (Int): Number of points 
    tol (Float64): tolerance of output (lower is more accurate)

    Outputs:
    X (Vector): Computed initial X values 
    Y (Vector): Computed initial Y values
    ϕ (Vector): Computed initial ϕ values
    =#
    zs,ws,PP = SSGW(kd,kH2,N,tol=tol)
    X = real(zs);
    Y = imag(zs);
    c = sqrt(9.81)*PP["cs"];
    ϕ = c*imag(hilbert(imag(zs)))
    return X,Y,ϕ,c    
end

function packet_manual(S,N,L,Δ,modes,k0,xbreaking,xinitial,left_lim,right_lim,sharpness; H=0)
    #=
    packet is a function that computes the initial X, Y, ϕ vectors for a focusing wave packet.
    This code is based off of Nick Pizzo's MATLAB code.

    Inputs:
    S (Float64): Linear prediction of maximum slope during focusing
    N (Int): Number of particles
    L (Float64): Length of periodic domain 
    Δ (Float64): Normalized bandwidth Δ = (k_max - k_0)/(2 k_0) where k_0 is central wavenumber
    modes (Int): Number of individual waves in packet 
    k0 (Float64): Central wavenumber

    Outputs:
    x_f (Vector): Output X
    y_f (Vector): Output Y
    f_f (Vector): Output ϕ
    =#

    ML = L; # Physical length of channel
    xₒ = 0;  
    x =collect( xₒ : ML/N : xₒ + (ML) * (1 - 1/N)); # domain
    wl = ML; # wl parameter is used in Dold code
    g = 9.81; # acc due to gravity. make sure Dold is set to same value
    k = k0; # central wavenumber
    w = sqrt(g*k); # associated central angular freq - deep water
    cg = g / 2 / w; # group velocity 
    D = Δ; # bandwidth (frequency)
    #xb = 20/D; #focusing location
    tb = (xbreaking-xinitial)/cg; # focusing time
    M = modes; # number of modes
    I = im; # define imaginary number
    kn = zeros(M); # allocate space for wavenumbers
    wn = zeros(M); # allocate space for frequencies
    an = zeros(M); # allocate space for amplitudes
    if ~iszero(H)
        w = sqrt(g*k*tanh(k*H))
        cg = 1/2 * w / k * (1 + 2*k*H/sinh(2*k*H))
        tb = (xbreaking-xinitial)/cg
    end
    for i ∈ 1 : M
        wn[i] = w *(1 - D *(i - M/2)/M)
        kn[i] = wn[i].^2 ./ g
        if ~iszero(H)
            kn[i] = invert_dispersion_relation(wn[i],H)
        end
        #kn[i] = k * ( 1 + D * (i - M/2) / M ); # define wavenumbers
        #wn[i] = sqrt( g * kn[i] ); # define associated freq
        an[i] = S / kn[i] / M; # define amplitudes
    end
    # Define surface displacement eta. Linear dispersive focusing 
    eta1 = zeros(length(x));
    for i ∈ 1 : M
        eta1 .+= an[i] * cos.( kn[i] .* (x .- (xbreaking)) .- wn[i] .* (0 - tb) );
    end

    # now window the initial waveform so we have 1 wave group

    # This will be updated, it's proportional to delta k. 
    #filt1 = 1/2 * (tanh.(0.28*(x.- -15))-tanh.(0.28.*(x.-10))); 
    filt1 = window_manual(x,left_lim,right_lim,sharpness,k)
    #filt1 = sech.(1/(2*Δ*k) .*(x.-xi));
    eta2 = eta1.*filt1;

    ## we now want to find the velocity potential
    ϕ = getϕ(x,eta2);
 
    x_f = x .- xₒ; # shift x-axis so it goes from 0 to wl    
    y_f = vec(eta2);
    f_f = ϕ;
    return x_f,y_f,f_f
end

function window_manual(a,left_lim,right_lim,sharpness,k₀)
    #=
    window is a function that computes the window as a function of time to isolate a single wave packet.

    =#
    

    return 0.5*(tanh.(sharpness/k₀ .* ((a .- left_lim))) .- tanh.(sharpness/k₀ .*( a .- right_lim)))


end

function getϕ(x,y;g=9.81)
    # This function finds ϕ from y using linear theory, FFTs, etc.

    # find spatial frequency 
    spatial_freq = 1/(x[2]-x[1]);
    # Compute fft of y (knowing it's real)
    transform = rfft(y);
    kvals = 2π.*rfftfreq(length(y),spatial_freq);

    # Multiply by relevant values 
    ϕtransform = transform[2:end].* (-im) .*(sqrt.(g./kvals[2:end]))
    # append zero to first value since zero mean 
    pushfirst!(ϕtransform,0);

    # Take inverse fft to recover ϕ
    return irfft(ϕtransform,length(y))

end