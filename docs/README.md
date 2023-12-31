# Castawave
![](breakingstepscropped.png)

Castawave is a boundary element solver for free surface potential flows, specifically adadpted to handling nonlinear surface waves. The mixed Euler-Lagrangian method used is inspired by [Dold (1992)](https://www.sciencedirect.com/science/article/pii/002199919290327U). 

This code base is using the [Julia Language](https://julialang.org/) and the
[DrWatson package](https://juliadynamics.github.io/DrWatson.jl/stable/) to make this scientific project and its results fast, easy to use, and reproducible.


It is authored by [Aidan Blaser](https://github.com/aidanblaser) & [Raphael Benamran](https://github.com/rbenamran).

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "Castawave"
```
which auto-activate the project and enable local path handling from DrWatson.
