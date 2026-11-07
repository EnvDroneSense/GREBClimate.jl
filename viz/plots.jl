# The four static plots and the helpers they share.

const VAR_INFO = Dict(
    :Ts => (label="Surface temperature", unit="K"),
    :Ta => (label="Air temperature", unit="K"),
    :To => (label="Deep-ocean temperature", unit="K"),
    :q => (label="Specific humidity", unit="kg/kg"),
    :albedo => (label="Albedo", unit=""),
    :ice => (label="Ice cover", unit="fraction"),
    :precip => (label="Precipitation", unit="mm/day"),
    :evap => (label="Evaporation", unit="mm/day"),
    :qcrcl => (label="Moisture circulation", unit="kg/kg"),
    :sw => (label="Net shortwave", unit="W/m²"),
    :lw => (label="Surface longwave", unit="W/m²"),
    :qlat => (label="Latent heat flux", unit="W/m²"),
    :qsens => (label="Sensible heat flux", unit="W/m²"),
    :olr => (label="Outgoing longwave", unit="W/m²"),
    :lwdown => (label="Downward longwave", unit="W/m²")
)

"Label and unit of a record field; unknown fields get their bare name."
fieldinfo(v::Symbol) = get(VAR_INFO, v, (label=String(v), unit=""))

unitlabel(v) = fieldinfo(v).unit == "" ? fieldinfo(v).label : "$(fieldinfo(v).label) [$(fieldinfo(v).unit)]"

"""
    runs(x) -> Vector of (; title, records, anomaly::Bool, times, dated::Bool)

The non-empty runs of a `greb_model!` result, or one untitled run for a bare
record vector. The scenario is an anomaly when the result's `scnr_anomaly` is
true; a result without that entry, and a bare record vector, count as absolute.

`times` is the `(year, month)` of each record: the result's `ctrl_time` and
`scnr_time`, and then `dated` is true. Records without times are counted from
January of year 1 ([`by_position`](@ref)).
"""
function runs(res::NamedTuple)
    out = [tagged(name, res[run], get(res, Symbol(run, :_time), nothing),
               run === :scnr && get(res, :scnr_anomaly, false))
           for (name, run) in (("control", :ctrl), ("scenario", :scnr)) if !isempty(res[run])]
    isempty(out) ? error("the result has no control or scenario records") : out
end
runs(records::AbstractVector) = [(title="", records=records, anomaly=false, times=by_position(length(records)), dated=false)]

function tagged(name, recs, times, anom)
    dated = times !== nothing
    return (title=anom ? "$name anomaly" : name, records=recs, anomaly=anom,
        times=dated ? times : by_position(length(recs)), dated=dated)
end

"Panels stacked in one column, each `h` pixels tall."
stack(panels, h) = length(panels) == 1 ? only(panels) :
                   Plots.plot(panels...; layout=(length(panels), 1), size=(760, h * length(panels)))

"Diverging colours centred on zero for anomalies, viridis otherwise."
colours(g, anom) = anom ? (c=:balance, clims=(-1, 1) .* max(maximum(abs, g), eps())) : (c=:viridis,)

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

"One line panel per run of `f(records, times)`, with a dashed zero line on anomalies."
function lines(f, x, var, xlabel; kw...)
    panels = map(runs(x)) do (t, recs, anom, times)
        y = f(recs, times)
        p = Plots.plot(1:length(y), y; xlabel=xlabel, ylabel=unitlabel(var), title=paneltitle(t, var),
            legend=false, lw=2, kw...)
        anom ? Plots.hline!(p, [0]; color=:gray, ls=:dash) : p
    end
    stack(panels, 280)
end

"""
    plot_timeseries(x; var=:Ts, annual=false, mask=nothing)

Area-weighted global mean of `var` per month, or per year with `annual=true`.
With a `mask` ([`land_mask`](@ref), [`band_mask`](@ref)) the mean is over its
cells.
"""
plot_timeseries(x; var::Symbol=:Ts, annual::Bool=false, mask=nothing) = lines(x, var, annual ? "year" : "month") do recs, times
    s = series(recs, var; mask=mask)
    annual ? GREBViz.annual(s, times) : s
end

"""
    plot_zonal(x; var=:Ts, month=:mean)

Zonal mean of `var` against latitude, one panel per run. `month` is as in
[`plot_map`](@ref).
"""
function plot_zonal(x; var::Symbol=:Ts, month=:mean)
    panels = map(runs(x)) do (t, recs, anom)
        z = zonal_mean(recs, var; month=month)
        p = Plots.plot(z, lats(length(z)); xlabel=unitlabel(var), ylabel="latitude [°N]",
            yticks=-90:30:90, title=paneltitle(t, var), legend=false, lw=2)
        anom ? Plots.vline!(p, [0]; color=:gray, ls=:dash) : p
    end
    stack(panels, 300)
end

"""
    plot_seasonal(x; var=:Ts)

The 12-month global-mean climatology, averaged over every whole year in the run.
"""
plot_seasonal(x; var::Symbol=:Ts) = lines((recs, times) -> seasonal_cycle(recs, var, times), x, var, "month"; marker=:circle, ms=3,
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
