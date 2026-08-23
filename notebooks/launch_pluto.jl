# Launch Pluto with the `viz` environment, which carries Plots, PlutoUI and a
# dev'd GREBClimate. Run from the repo root:
#
#   julia notebooks/launch_pluto.jl
#
# then open notebooks/GREB_explorer.jl in the browser tab it opens.
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "viz"))
Pkg.instantiate()

using Pluto
Pluto.run(notebook = joinpath(@__DIR__, "GREB_explorer.jl"))
