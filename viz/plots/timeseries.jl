# Evolution over time. Also the template for any new visualization: `compute`
# returns numbers, `plot` returns a figure, `register!` wires it in.
#
# Both plots here carry each run's `kind` through to the plot, so absolutes and
# anomalies land on separate panels rather than one axis where a 2 K signal
# disappears against a 288 K control.

register!(
    key   = :timeseries,
    label = "Evolution over time",
    kind  = :series,
    help  = "Area-weighted regional mean of one variable, month by month, for " *
            "every selected run. Absolute and anomaly runs get their own panel.",

    compute = (runs, o) -> map(runs) do r
        s = series(records_of(r), o.var; region=o.region, fields=o.fields)
        (label_of(r), kind_of(r), o.annual ? annual(s) : s)
    end,

    plot = function (data, o)
        info = varinfo(o.var)
        xlab = o.annual ? "year" : "month"
        panels = map(by_kind(data)) do (k, group)
            p = Plots.plot(; xlabel=xlab, ylabel=unitlabel(info) * kind_suffix(k),
                           title="$(info.label)$(kind_suffix(k)) - $(regionlabel(o.region))",
                           legend=:outertopright)
            for (name, y) in group
                isempty(y) && continue
                Plots.plot!(p, 1:length(y), y; label=name, lw=2)
            end
            k === :anomaly && Plots.hline!(p, [0]; color=:gray, ls=:dash, label="")
            p
        end
        length(panels) == 1 ? panels[1] :
            Plots.plot(panels...; layout=(length(panels), 1), size=(760, 280 * length(panels)))
    end,
)

register!(
    key   = :seasonal,
    label = "Seasonal cycle",
    kind  = :series,
    help  = "The 12-month climatology, averaged over every whole year in the run.",

    compute = (runs, o) -> map(runs) do r
        (label_of(r), kind_of(r),
         seasonal_cycle(records_of(r), o.var; region=o.region, fields=o.fields))
    end,

    plot = function (data, o)
        info = varinfo(o.var)
        panels = map(by_kind(data)) do (k, group)
            p = Plots.plot(; xlabel="month", ylabel=unitlabel(info) * kind_suffix(k),
                           title="Seasonal cycle$(kind_suffix(k)) - $(info.label), $(regionlabel(o.region))",
                           xticks=(1:12, ["J","F","M","A","M","J","J","A","S","O","N","D"]),
                           legend=:outertopright)
            for (name, y) in group
                isempty(y) && continue
                Plots.plot!(p, 1:12, y; label=name, lw=2, marker=:circle, ms=3)
            end
            p
        end
        length(panels) == 1 ? panels[1] :
            Plots.plot(panels...; layout=(length(panels), 1), size=(760, 280 * length(panels)))
    end,
)
