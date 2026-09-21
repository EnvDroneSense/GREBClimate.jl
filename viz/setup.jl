# One-time setup of the viz/ environment: link it to this clone of GREBClimate
# (not registered) and install Plots, Pluto and PlutoUI.
#
#   julia viz/setup.jl

using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.instantiate()
