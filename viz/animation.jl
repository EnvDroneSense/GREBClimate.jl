# Maps through a run: one frame per month or year, runs side by side, each on a
# colour scale fixed across its frames so a colour change is a value change.

const MONTH_ABBR = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

"""
    frame_limits(frames, anom) -> (lo, hi)

Colour limits shared by every frame. Absolute: the full range. Anomaly: symmetric
about zero at the 99th percentile of |value|, because a few cells reach ±23 K
where most sit at 1-3 K (15-year 2xCO2 run) and would wash the rest out.
"""
function frame_limits(frames, anom::Bool)
    if anom
        m = quantile(abs.(reduce(vcat, vec.(frames))), 0.99)
        return m == 0 ? (-1.0, 1.0) : (-m, m)
    end
    lo, hi = minimum(minimum, frames), maximum(maximum, frames)
    lo == hi ? (lo - 1, hi + 1) : (lo, hi)
end

"Frames and colour limits per run, computed once so drawing a frame only draws."
struct Evolution
    var::Symbol
    step::Symbol
    panels::Vector{NamedTuple}
end

"""
    evolution(x; var=:Ts, step=:month) -> Evolution

Frames of `var` for each run of a `greb_model!` result (or a record vector), per
month or per whole year (`step=:year`).
"""
function evolution(x; var::Symbol=:Ts, step::Symbol=:month)
    panels = map(runs(x)) do (t, recs, anom)
        f = map_frames(recs, var; step=step)
        isempty(f) && error("run '$t' has no whole $step to show")
        (title=t, frames=f, anom=anom, clims=Float64.(frame_limits(f, anom)))
    end
    Evolution(var, step, panels)
end

"Frames in the longest run."
frame_count(ev::Evolution) = maximum(length(p.frames) for p in ev.panels)

"""
    evolution_frame(ev, i; fields=nothing)

Frame `i` of every run, side by side. A shorter run cycles (`mod1(i, n)`), so
with monthly steps all panels show the same calendar month.
"""
function evolution_frame(ev::Evolution, i::Integer; fields=nothing)
    panels = map(ev.panels) do p
        n = length(p.frames)
        k = mod1(i, n)
        when = ev.step === :year ? "year $k" : "$(MONTH_ABBR[mod1(k, 12)]), year $(cld(k, 12))"
        geomap(p.frames[k]; fields=fields, c=p.anom ? :balance : :viridis, clims=p.clims,
               title=lstrip("$(p.title): $when  ($k/$n)", [':', ' ']), titlefontsize=10,
               colorbar_title=unitlabel(ev.var))
    end
    length(panels) == 1 ? only(panels) :
        Plots.plot(panels...; layout=(1, length(panels)), size=(600 * length(panels), 340),
                   left_margin=7Plots.mm, bottom_margin=5Plots.mm)
end

"""
    evolution_gif(path, ev; fps=6, fields=nothing) -> Plots.AnimatedGif

Write every frame to a GIF at `path`; the result displays inline in Pluto or Jupyter.
"""
function evolution_gif(path::AbstractString, ev::Evolution; fps::Real=6, fields=nothing)
    anim = Plots.Animation()
    foreach(i -> Plots.frame(anim, evolution_frame(ev, i; fields=fields)), 1:frame_count(ev))
    Plots.gif(anim, path; fps=fps, show_msg=false)
end
