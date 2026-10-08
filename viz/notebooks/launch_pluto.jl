# Open the GREB explorer notebook in Pluto. Run from the repository root:
#
#   julia viz/notebooks/launch_pluto.jl

using Pkg
viz = joinpath(@__DIR__, "..")
isfile(joinpath(viz, "Manifest.toml")) || include(joinpath(viz, "setup.jl"))   # first launch only
Pkg.activate(viz)

using Pluto
Pluto.run(notebook = joinpath(@__DIR__, "GREB_explorer.jl"))
