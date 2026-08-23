# Spatial maps, on real lon/lat axes rather than raw grid indices.

register!(
    key   = :map,
    label = "Map (spatial field)",
    kind  = :map,
    help  = "One variable over the globe. `month=:mean` averages the whole run, " *
            "which is usually what you want for anything longer than a year.",

    compute = (runs, o) -> map(runs) do r
        (label_of(r), kind_of(r), field(records_of(r), o.var; month=o.month))
    end,

    plot = function (data, o)
        info = varinfo(o.var)
        panels = map(data) do (name, k, g)
            nx, ny = size(g)
            kindtag = kind_suffix(k)
            pnl = Plots.heatmap(lons(nx), lats(ny), permutedims(g);
                xlabel="longitude [°E]", ylabel="latitude [°N]",
                title="$name - $(info.label)$(kindtag)", c=:viridis,
                colorbar_title=info.unit, yticks=-90:30:90, xticks=0:60:360)
            o.coastlines && coastlines!(pnl, o.fields)
            pnl
        end
        length(panels) == 1 ? panels[1] :
            Plots.plot(panels...; layout=(length(panels), 1),
                       size=(760, 300 * length(panels)))
    end,
)

register!(
    key   = :diffmap,
    label = "Difference map (A - B)",
    kind  = :diffmap,
    help  = "Second run subtracted from the first, on a diverging scale centred " *
            "at zero. This is how you see WHERE a scenario acts, not just by how much.",

    compute = function (runs, o)
        length(runs) < 2 && error("difference map needs two runs; got $(length(runs))")
        a, b = runs[1], runs[2]
        kind_of(a) === kind_of(b) || error("""
            cannot difference a :$(kind_of(a)) run against a :$(kind_of(b)) one -
            '$(label_of(a))' and '$(label_of(b))' are not the same kind of quantity.
            A scenario from greb_model! is already an anomaly against its control,
            so subtracting the control from it is meaningless. Either compare two
            runs of the same kind, or convert first with
            `to_absolute(scenario_run, control_run)`.""")
        (label_of(a), label_of(b), kind_of(a),
         difference(records_of(a), records_of(b), o.var; month=o.month))
    end,

    plot = function (data, o)
        aname, bname, k, d = data
        info = varinfo(o.var)
        nx, ny = size(d)
        lim = maximum(abs, d)
        lim = lim == 0 ? 1.0 : lim
        p = Plots.heatmap(lons(nx), lats(ny), permutedims(d);
            xlabel="longitude [°E]", ylabel="latitude [°N]",
            title="$aname - $bname  ($(info.label)$(kind_suffix(k)))",
            c=:balance, clims=(-lim, lim), colorbar_title=info.unit,
            yticks=-90:30:90, xticks=0:60:360)
        o.coastlines && coastlines!(p, o.fields)
        p
    end,
)
