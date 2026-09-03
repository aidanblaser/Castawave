# Castawave
![](breakingstepscropped.png)

Castawave is a boundary element solver for free surface potential flows, specifically adadpted to handling nonlinear surface waves. The mixed Euler-Lagrangian method used is inspired by [Dold (1992)](https://www.sciencedirect.com/science/article/pii/002199919290327U).

This code base uses the [Julia Language](https://julialang.org/). It is a
standard, installable Julia package (it has its own `Project.toml` with a
UUID and version) rather than a DrWatson-style research project - that
keeps it lightweight and easy to depend on from other projects, including
your own analysis scripts.

It is authored by [Aidan Blaser](https://github.com/aidanblaser) & [Raphael Benamran](https://github.com/rbenamran).

To use it locally, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.develop(path="path/to/this/project")
   julia> using Castawave
   ```
   `Pkg.develop` registers Castawave as a local dependency of whatever
   project/environment is currently active, without copying or vendoring
   it - edits you make to this repo are picked up immediately.

See `scripts/QuickStart.jl` for a minimal, self-contained worked example
(it ships its own nested environment, `scripts/Project.toml`, so it needs
no manual setup beyond the one-time `Pkg.develop` + `Pkg.instantiate` shown
in that script's header comment).

