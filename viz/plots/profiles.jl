# Latitude-resolved views: the conventional way climate output gets read.

register!(
    key   = :zonal,
    label = "Zonal-mean profile",
    kind  = :profile,
    help  = "Longitude averaged away, leaving value against latitude. The " *
            "quickest way to see polar amplification or a tropical bias.",

    compute = (runs, o) -> map(runs) do r
        (label_of(r), kind_of(r), zonal_mean(records_of(r), o.var; month=o.month))
    end,

    plot = function (data, o)
        info = varinfo(o.var)
        panels = map(by_kind(data)) do (k, group)
            p = Plots.plot(; xlabel=unitlabel(info) * kind_suffix(k), ylabel="latitude [°N]",
                           title="Zonal mean$(kind_suffix(k)) - $(info.label)",
                           yticks=-90:30:90, legend=:outertopright)
            for (name, v) in group
                isempty(v) && continue
                Plots.plot!(p, v, lats(length(v)); label=name, lw=2)
            end
            k === :anomaly && Plots.vline!(p, [0]; color=:gray, ls=:dash, label="")
            p
        end
        length(panels) == 1 ? panels[1] :
            Plots.plot(panels...; layout=(1, length(panels)), size=(420 * length(panels), 420))
    end,
)

register!(
    key   = :hovmoller,
    label = "Hovmöller (time × latitude)",
    kind  = :hovmoller,
    help  = "Zonal mean plotted against both time and latitude, so you can watch " *
            "a signal migrate. Shows the first selected run only.",

    compute = function (runs, o)
        isempty(runs) && error("Hovmöller needs at least one run")
        r = runs[1]
        t, φ, m = hovmoller(records_of(r), o.var)
        (label_of(r), t, φ, m)
    end,

    plot = function (data, o)
        name, t, φ, m = data
        info = varinfo(o.var)
        Plots.heatmap(collect(t), φ, permutedims(m);
            xlabel="month", ylabel="latitude [°N]",
            title="$name - $(info.label)", c=:viridis,
            colorbar_title=info.unit, yticks=-90:30:90)
    end,
)
