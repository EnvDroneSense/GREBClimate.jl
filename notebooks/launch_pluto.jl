# Open the GREB explorer notebook in Pluto. Run from the repository root:
#
#   julia notebooks/launch_pluto.jl

using Pkg
viz = joinpath(@__DIR__, "..", "viz")
isfile(joinpath(viz, "Manifest.toml")) || include(joinpath(viz, "setup.jl"))   # first launch only
Pkg.activate(viz)

using Pluto
Pluto.run(notebook = joinpath(@__DIR__, "GREB_explorer.jl"))
