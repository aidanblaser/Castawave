#=
Castawave.jl is the module entry point for the Castawave package: a
boundary-element solver for nonlinear free-surface water waves, adapted
from Dold (1992, J. Comp. Phys. 103, 90-115).

Load it with:

    using DrWatson
    @quickactivate "Castawave"
    include(projectdir()*"/src/Castawave.jl")
    using .Castawave

See scripts/QuickStart.jl for a minimal worked example of exactly that.

ClamondIC.jl and DoldFcns.jl (an initial-condition generator and a
Dold-I/O-format-compatible validation driver, respectively) have been
removed from this package - they were specific to other, unrelated
analysis rather than the core solver. They're still recoverable from git
history if ever needed again.
=#

module Castawave

using DrWatson
@quickactivate "Castawave"

# MainSolver.jl's own header already includes Constants.jl, Types.jl and
# HelperFunctions.jl (in that order - Types.jl has to precede
# HelperFunctions.jl; see the note at the top of Types.jl) and declares
# its own package dependencies. Reusing that existing chain here, rather
# than duplicating it, means this file can't drift out of sync with it.
include("MainSolver.jl")

export SimulationParameters, SolverWorkspace, runSim,
       fixedTimeOperations, fixedTimeOperations!

end # module Castawave
