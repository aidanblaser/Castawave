using DrWatson, Test
@quickactivate "Castawave"

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "Castawave tests" begin
    @test 1 == 1
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")







function DDI1(Ω::Vector, N, q=0)
    #=
    DDI1 is a function that implements the 11-point Lagrangian polynomial interpolation for the first coefficient from eq... , therefore returning the first order derivative value at the center point.
    =#
    if q == 0
        Ω_p = zeros(Complex, N)
    elseif q == 1
        Ω_p = zeros(N)
    end

    Threads.@threads for i in 1:N
        points = [Ω[mod1(i+1, N)] - Ω[mod1(i-1, N)], Ω[mod1(i+2, N)] - Ω[mod1(i-2, N)], Ω[mod1(i+3, N)] - Ω[mod1(i-3, N)], Ω[mod1(i+4, N)] - Ω[mod1(i-4, N)], Ω[mod1(i+5, N)] - Ω[mod1(i-5, N)]]
        Ω_p[i] = sum([1050.0, -300.0, 75.0, -12.5, 1.0] ./ 1260.0 .* points)
    end

    return Ω_p
end


L = 128
x = [L / 2 / π * sin(2 * π * i / L) for i in 0:L-1]
@time y = DDI1(x, L)
y[1]
