#= 
This file defines constants for the program that are .
=#
using JLD2
using Plots

# df = jldopen("0.01Doldsim.jld2")
# keys(df)

# X = df["x_o"]
# Y = df["y_o"]
# T = df["t"]

# c = union(1, collect(3:22))
# function fn()
#     anim = @animate for i in 1:101
#         scatter([X[1,:]], [Y[1,:]], label = "Timestepped", legend = :topright, framestyle= :box,background_color="black", markerstrokewidth=0,markersize=1, dpi = 300, xlabel=L"x \,(m)",ylabel=L"\eta \,(m)", title= @sprintf("Time: %.2f s", T[i]))
#         scatter!([X[i,:]], [Y[i,:]], label = "Initial position", framestyle= :box,background_color="black", markerstrokewidth=0,markersize=1)
#     end every 1
#     gif(anim, "doldlinear.gif", fps=4)
# end
# fn()

FACTOR1 = [2100 -600 150 -25 2; 
-70098 52428 -14607 2522 -205;
1938 -1872 783 -152 13;
-378 408 -207 52 -5;
42 -48 27 -8 1]

WEIGHTS1 = diagm([2520, 181440, 34560, 120960, 725760])

FACTOR2 = [-73766 42000 -6000 1000 -125 8;
192654 -140196 52428 -9738 1261 -82;
-12276 9690 -4680 1305 -190 13;
462 -378 204 -69 13 -1;
-252 210 -120 45 -10 1]

WEIGHTS2 = diagm([50400, 362880, 172800, 120960, 3628800])

GRAVITY = 9.81

COEFFICIENTS1 = [1050.0, -300.0, 75.0, -12.5, 1.0] ./ 1260.0

COEFFICIENTS2 = [-9220.75, 5250.0, -750.0, 125.0, -15.625, 1.0] ./ 3150.0

SMOOTHCOEFFICIENTS = [3432, -3003, 2002, -1001, 364, -91, 14, -1]
