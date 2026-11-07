# Plots that put runs side by side: the difference of two, several as lines, and
# one number per ensemble member.

"""
    plot_difference(a, b; var=:Ts, month=:mean, fields=nothing)

`a - b` for two record vectors of the same kind (both absolute or both
anomalies), on a scale centred at zero: where one run differs from the other,
not just by how much. For a scenario against its control:
`plot_difference(to_absolute(res.scnr, res.ctrl), res.ctrl)`; for two
scenarios of the same preset: `plot_difference(r1.scnr, r2.scnr)`.
"""
function plot_difference(a::AbstractVector, b::AbstractVector; var::Symbol=:Ts, month=:mean, fields=nothing)
    d = field(a, var; month=month) .- field(b, var; month=month)
    geomap(d; fields=fields, title="Difference - $(fieldinfo(var).label)",
        colorbar_title=fieldinfo(var).unit, colours(d, true)...)
end

"""
    plot_compare(members; var=:Ts, annual=true, mask=nothing)

The global mean of `var` for several runs on one set of axes. `members` is a
vector of `label => result` pairs; each result is a `greb_model!` result (its
scenario is drawn, or its control when it has no scenario) or a record vector.
Anomalies and absolute runs should not be mixed.
"""
function plot_compare(members::AbstractVector{<:Pair}; var::Symbol=:Ts, annual::Bool=true, mask=nothing)
    isempty(members) && error("no members to compare")
    p = Plots.plot(; xlabel=annual ? "year" : "month", ylabel=unitlabel(var), lw=2, legend=:outertopright)
    zero_line = false
    for (label, x) in members
        run = last(runs(x))
        s = series(run.records, var; mask=mask)
        y = annual ? GREBViz.annual(s, run.times) : s
        Plots.plot!(p, 1:length(y), y; label=string(label))
        zero_line |= run.anomaly
    end
    zero_line ? Plots.hline!(p, [0]; color=:gray, ls=:dash, label="") : p
end

"""
    plot_members(members; xlabel="", baseline=nothing)

One bar per ensemble member. `members` is a vector of `label => number` pairs,
for example the warming in the last year of each `run_ensemble` member. With
`baseline` set to one of the labels, the bars show the difference to that
member.
"""
function plot_members(members::AbstractVector{<:Pair}; xlabel::AbstractString="", baseline=nothing)
    isempty(members) && error("no members to plot")
    labels = [string(first(m)) for m in members]
    values = Float64[last(m) for m in members]
    if baseline !== nothing
        k = findfirst(==(string(baseline)), labels)
        k === nothing && error("baseline '$baseline' is not among the members: $(join(labels, ", "))")
        values .-= values[k]
        xlabel = isempty(xlabel) ? "difference to $baseline" : xlabel
    end
    n = length(labels)
    lo, hi = min(0.0, minimum(values)), max(0.0, maximum(values))
    pad = (hi - lo == 0 ? 1.0 : hi - lo) * 0.08
    Plots.bar(1:n, values; orientation=:h, legend=false, xlabel=xlabel, yticks=(1:n, labels),
        yflip=true, ylims=(0.4, n + 0.6), xlims=(lo - pad, hi + pad),
        size=(760, 110 + 36 * n), left_margin=8Plots.mm, bottom_margin=5Plots.mm)
end
