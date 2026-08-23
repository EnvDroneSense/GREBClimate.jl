"""
    GREBViz

Visualization tooling for `greb_model!` output. **Not part of the GREBClimate
package** - a separate tool that reads what the model returns, so it adds no
dependency to the package itself.

    include("viz/GREBViz.jl"); using .GREBViz

    res = greb_model!(RunSpec(ctrl=10, scnr=30), cfg; jld2_dir=dir, fields=fields)
    runs = ["control" => res.ctrl, "2xCO2" => res.scnr]

    render(:timeseries, runs; var=:Ts, region=:arctic, annual=true)
    render(:diffmap,    runs; var=:Ts)

# Adding a visualization

Drop a file in `viz/plots/` that calls [`register!`](@ref). It is picked up
automatically - neither this file nor the notebook needs editing. See
`viz/plots/timeseries.jl` for the shortest example.

Two layers, kept apart on purpose:

  * **Reductions** (`reductions.jl`) - pure, no plotting, usable from any
    script and testable on their own.
  * **Visualizations** (`plots/`) - each declares `compute` (numbers) and
    `plot` (picture) separately, so styling and maths change independently.

Regional means are **area-weighted by `cos(latitude)`** by default; see
[`area_weights`](@ref) for why that matters.
"""
module GREBViz

import Plots

export render, numbers, vizlist, viz_options, var_options, region_options, variables,
       series, annual, zonal_mean, seasonal_cycle, hovmoller, field, difference,
       area_weights, region_mask, lats, lons, register!, register_region!,
       varinfo, VizSpec, Run, control, scenario, greb_runs, to_absolute,
       monthly_climatology, mixed_kinds, coastlines!

include("registry.jl")
include("runs.jl")
include("reductions.jl")

# Auto-include every visualization. Sorted so load order is reproducible.
for _f in sort(readdir(joinpath(@__DIR__, "plots"); join=true))
    endswith(_f, ".jl") && include(_f)
end

# ── Presentation helpers shared by the plot files ──────────────────────────
unitlabel(info) = info.unit == "" ? info.label : "$(info.label) [$(info.unit)]"

"""
    coastlines!(p, fields; color=:black, lw=0.7)

Draw the land outline over an existing map by contouring `fields.z_topo` at the
0 m shoreline. No-op when `fields` is `nothing`, so a map without model data
still renders.

Only `.z_topo` is touched, so anything carrying that field works - the tool
stays independent of GREBClimate's types.
"""
function coastlines!(p, fields; color=:black, lw=0.7)
    fields === nothing && return p
    zt = fields.z_topo
    nx, ny = size(zt)
    Plots.contour!(p, lons(nx), lats(ny), permutedims(Float64.(zt));
                   levels=[0.0], linecolor=color, linewidth=lw,
                   colorbar_entry=false, label="")
    p
end

function regionlabel(key::Symbol)
    i = findfirst(p -> p.first === key, REGIONS)
    i === nothing ? String(key) : REGIONS[i].second.label
end

# ── Entry point ────────────────────────────────────────────────────────────

"""
    default_opts(; kw...)

Every option any visualization understands, with defaults. Passing an option a
given plot ignores is harmless, which keeps the notebook's control panel
uniform across plot types.
"""
default_opts(; var::Symbol=:Ts, region::Symbol=:global, fields=nothing,
               annual::Bool=false, month=:mean, coastlines::Bool=true) =
    (var=var, region=region, fields=fields, annual=annual, month=month,
     coastlines=coastlines)

"""
    normalize_runs(runs)

Accepts a `Run`, `res.ctrl`, `"name" => res.ctrl`, or a vector of any of those,
and always returns `Vector{Run}`. Untagged input is assumed `:absolute` - use
[`greb_runs`](@ref) to tag a `greb_model!` result correctly in one call.
"""
function normalize_runs(runs)
    runs isa Run && return [runs]
    runs isa Pair && return [Run(runs.first, runs.second)]
    isempty(runs) && return Run[]
    f = first(runs)
    (f isa Run || f isa Pair) || return [Run("run", runs)]
    [x isa Run ? x : Run(label_of(x), records_of(x)) for x in runs]
end

"""
    by_kind(data)

Group `(label, kind, values)` triples by kind, preserving order. Plots use this
to put absolutes and anomalies on separate panels instead of one useless axis.
"""
function by_kind(data)
    ks = unique(k for (_, k, _) in data)
    [(k, [(l, v) for (l, kk, v) in data if kk === k]) for k in ks]
end

"""
    render(key, runs; var=:Ts, region=:global, annual=false, month=:mean, fields=nothing)

Compute and plot visualization `key`. `runs` is one run or several, optionally
labelled: `["control" => res.ctrl, "2xCO2" => res.scnr]`.

`month` accepts an index, `:last`, or `:mean` (time-mean over the run).
`fields` is needed for the `:land`/`:ocean` regions and for the coastline
overlay on maps, both of which read `z_topo`; pass `coastlines=false` to skip
the outline.
"""
function render(key::Symbol, runs; kw...)
    spec = viz(key)
    o = merge(default_opts(), NamedTuple(kw))
    rs = normalize_runs(runs)
    isempty(rs) && error("no runs given")
    vars = variables(records_of(rs[1]))
    o.var in vars || error("unknown variable :$(o.var); records carry $(join(vars, ", "))")
    spec.plot(spec.compute(rs, o), o)
end

"""
    numbers(key, runs; kw...)

The `compute` half of a visualization without plotting it - for when you want
the values rather than the picture.
"""
numbers(key::Symbol, runs; kw...) =
    viz(key).compute(normalize_runs(runs), merge(default_opts(), NamedTuple(kw)))

end # module
