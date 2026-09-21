# The four static plots and the helpers they share.

const VAR_INFO = Dict(
    :Ts     => (label="Surface temperature",    unit="K"),
    :Ta     => (label="Air temperature",        unit="K"),
    :To     => (label="Deep-ocean temperature", unit="K"),
    :q      => (label="Specific humidity",      unit="kg/kg"),
    :albedo => (label="Albedo",                 unit=""),
    :ice    => (label="Ice cover",              unit="fraction"),
    :precip => (label="Precipitation",          unit="mm/day"),
    :evap   => (label="Evaporation",            unit="mm/day"),
    :qcrcl  => (label="Moisture circulation",   unit="kg/kg"),
    :sw     => (label="Net shortwave",          unit="W/m²"),
    :lw     => (label="Surface longwave",       unit="W/m²"),
    :qlat   => (label="Latent heat flux",       unit="W/m²"),
    :qsens  => (label="Sensible heat flux",     unit="W/m²"),
)

"Label and unit of a record field; unknown fields get their bare name."
fieldinfo(v::Symbol) = get(VAR_INFO, v, (label=String(v), unit=""))

unitlabel(v) = fieldinfo(v).unit == "" ? fieldinfo(v).label : "$(fieldinfo(v).label) [$(fieldinfo(v).unit)]"

"""
    runs(x) -> Vector of (title, records, anomaly::Bool)

The non-empty runs of a `greb_model!` result, or one untitled run for a bare
record vector. `res.scnr` is usually an anomaly against the control's final
year, but stays absolute for orbital experiments and when `ctrl=0`, so the kind
is read from the data: absolute Ts averages ~280 K, an anomaly a few K.
"""
function runs(res::NamedTuple)
    out = [tagged(name, recs) for (name, recs) in (("control", res.ctrl), ("scenario", res.scnr))
           if !isempty(recs)]
    isempty(out) ? error("the result has no control or scenario records") : out
end
runs(records::AbstractVector) = [("", records, false)]

tagged(name, recs) = (anom = mean(first(recs).Ts) < 100; (anom ? "$name anomaly" : name, recs, anom))

"Panels stacked in one column, each `h` pixels tall."
stack(panels, h) = length(panels) == 1 ? only(panels) :
    Plots.plot(panels...; layout=(length(panels), 1), size=(760, h * length(panels)))

"Diverging colours centred on zero for anomalies, viridis otherwise."
colours(g, anom) = anom ? (c=:balance, clims=(-1, 1) .* max(maximum(abs, g), eps())) : (c=:viridis,)

"""
    coastlines!(p, fields; color=:black, lw=0.7)

Outline the cells where `fields.z_topo > 0`. No-op when `fields` is `nothing`.
Drawn as line segments, not a contour: a contour joins the heatmap's colour
scale and is clipped when fixed `clims` exclude 0 m.
"""
function coastlines!(p, fields; color=:black, lw=0.7)
    fields === nothing && return p
    land = fields.z_topo .> 0
    nx, ny = size(land)
    dx, dy = 360 / nx, 180 / ny
    xs, ys = Float64[], Float64[]
    seg!(x1, y1, x2, y2) = (append!(xs, (x1, x2, NaN)); append!(ys, (y1, y2, NaN)))
    for j in 1:ny, i in 1:nx
        y0, y1 = -90 + (j - 1) * dy, -90 + j * dy
        land[i, j] != land[mod1(i + 1, nx), j] && seg!(i * dx, y0, i * dx, y1)       # east edge, wraps
        j < ny && land[i, j] != land[i, j + 1] && seg!((i - 1) * dx, y1, i * dx, y1) # north edge
    end
    Plots.plot!(p, xs, ys; linecolor=color, linewidth=lw, label="", xlims=(0, 360), ylims=(-90, 90))
end

"Heatmap of a lon x lat grid on degree axes, outlined when `fields` is given."
function geomap(g::AbstractMatrix; fields=nothing, kw...)
    nx, ny = size(g)
    p = Plots.heatmap(lons(nx), lats(ny), permutedims(g); xlabel="longitude [°E]",
                      ylabel="latitude [°N]", xticks=0:60:360, yticks=-90:30:90, kw...)
    coastlines!(p, fields)
end

paneltitle(t, var) = isempty(t) ? fieldinfo(var).label : "$t - $(fieldinfo(var).label)"

"""
    plot_map(x; var=:Ts, month=:mean, fields=nothing)

`var` over the globe for each run: the time-mean (`:mean`), the last month
(`:last`) or a month index. Pass the model's `fields` to draw coastlines.
"""
function plot_map(x; var::Symbol=:Ts, month=:mean, fields=nothing)
    panels = map(runs(x)) do (t, recs, anom)
        g = field(recs, var; month=month)
        geomap(g; fields=fields, title=paneltitle(t, var), colorbar_title=fieldinfo(var).unit,
               colours(g, anom)...)
    end
    stack(panels, 300)
end

"One line panel per run of `f(records)`, with a dashed zero line on anomalies."
function lines(f, x, var, xlabel; kw...)
    panels = map(runs(x)) do (t, recs, anom)
        y = f(recs)
        p = Plots.plot(1:length(y), y; xlabel=xlabel, ylabel=unitlabel(var), title=paneltitle(t, var),
                       legend=false, lw=2, kw...)
        anom ? Plots.hline!(p, [0]; color=:gray, ls=:dash) : p
    end
    stack(panels, 280)
end

"""
    plot_timeseries(x; var=:Ts, annual=false)

Area-weighted global mean of `var` per month, or per year with `annual=true`.
"""
plot_timeseries(x; var::Symbol=:Ts, annual::Bool=false) =
    lines(x, var, annual ? "year" : "month") do recs
        s = series(recs, var)
        annual ? GREBViz.annual(s) : s
    end

"""
    plot_seasonal(x; var=:Ts)

The 12-month global-mean climatology, averaged over every whole year in the run.
"""
plot_seasonal(x; var::Symbol=:Ts) =
    lines(recs -> seasonal_cycle(recs, var), x, var, "month"; marker=:circle, ms=3,
          xticks=(1:12, ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]))

"""
    plot_hovmoller(x; var=:Ts)

Zonal mean of `var` against month and latitude, one panel per run.
"""
function plot_hovmoller(x; var::Symbol=:Ts)
    panels = map(runs(x)) do (t, recs, anom)
        months, φ, m = hovmoller(recs, var)
        Plots.heatmap(collect(months), φ, permutedims(m); xlabel="month", ylabel="latitude [°N]",
                      yticks=-90:30:90, colorbar_title=fieldinfo(var).unit, title=paneltitle(t, var),
                      colours(m, anom)...)
    end
    stack(panels, 300)
end
