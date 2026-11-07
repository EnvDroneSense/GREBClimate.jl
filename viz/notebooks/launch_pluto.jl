# Open the GREB explorer notebook in Pluto. Run from the repository root:
#
#   julia viz/notebooks/launch_pluto.jl

using Pkg
viz = joinpath(@__DIR__, "..")
include(joinpath(viz, "setup.jl"))   # seconds when the environment is current
Pkg.activate(viz)

using Pluto
Pluto.run(notebook=joinpath(@__DIR__, "GREB_explorer.jl"))
