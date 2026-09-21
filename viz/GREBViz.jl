"""
    GREBViz

Basic plots of `greb_model!` output. Not part of the GREBClimate package: it has
its own environment (`viz/Project.toml`).

    include("viz/GREBViz.jl"); using .GREBViz

    res = greb_model!(RunSpec(ctrl=5, scnr=15), cfg; jld2_dir=dir, fields=fields)
    plot_map(res; var=:Ts, fields=fields)
    plot_timeseries(res; annual=true)

Every plot takes a `greb_model!` result, drawing control and scenario as separate
panels, or a bare vector of monthly records.
"""
module GREBViz

import Plots
using Statistics: mean, quantile

export plot_map, plot_timeseries, plot_seasonal, plot_hovmoller, coastlines!,
       evolution, evolution_frame, evolution_gif, frame_count,
       series, annual, seasonal_cycle, field, hovmoller, map_frames,
       lats, lons, area_weights, fieldinfo

include("reductions.jl")
include("plots.jl")
include("animation.jl")

end # module
