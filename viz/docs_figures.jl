# Regenerate the figures of the docs pages (docs/src/assets/viz/*.png): runs a
# 2xCO2 scenario and saves the map, Hovmoller and difference plots. Needs the
# dataset; about a minute.
#
#   julia --project=viz -t 2,0 viz/docs_figures.jl [output_dir]    # default: docs/src/assets/viz

using Plots, GREBClimate
gr()
include(joinpath(@__DIR__, "GREBViz.jl"))
using .GREBViz

outdir = mkpath(get(ARGS, 1, joinpath(@__DIR__, "..", "docs", "src", "assets", "viz")))
dir = greb_data_dir(allow_download=false)
fields = load_climatology(dir; dataset=:ncep)

@info "running GREB: 2xCO2 (3 yr control, 15 yr scenario)"
# Base.CoreLogging: the viz environment does not list the Logging stdlib.
res = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
    greb_model!(RunSpec(ctrl=3, scnr=15), preset(:co2_double); jld2_dir=dir, fields=fields)
end

save(p, name) = (savefig(p, joinpath(outdir, "$name.png")); println("  -> $name.png"))
save(plot_map(res; fields=fields), "map")
save(plot_hovmoller(res), "hovmoller")
save(plot_difference(to_absolute(res.scnr, res.ctrl), res.ctrl; fields=fields), "difference")
