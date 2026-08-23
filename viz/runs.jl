# What a "run" is, and why it needs a kind.
#
# `greb_model!` returns two different kinds of quantity. `res.ctrl` is
# absolute (global-mean Ts ~288 K). `res.scnr` has already been through
# `apply_scenario_anomalies`, so it is an anomaly against the control
# climatology (~2 K). Nothing in the returned arrays distinguishes them, so
# plotting both on one axis silently flattens the entire scenario signal
# against a 288 K control - and a difference map computes 288 - 2 and saturates.

"""
    Run(label, records; kind=:absolute)

A labelled set of monthly records. `kind` is `:absolute` (real units, e.g. K)
or `:anomaly` (a difference against a control climatology).

Prefer [`control`](@ref)/[`scenario`](@ref), or [`greb_runs`](@ref) to build
both straight from a `greb_model!` result.
"""
struct Run
    label::String
    records::Vector
    kind::Symbol
end

function Run(label, records; kind::Symbol=:absolute)
    kind in (:absolute, :anomaly) ||
        error("kind must be :absolute or :anomaly, got :$kind")
    Run(String(label), collect(records), kind)
end

"A run in real units."
control(label, records) = Run(label, records; kind=:absolute)

"A run of anomalies against a control climatology."
scenario(label, records) = Run(label, records; kind=:anomaly)

records_of(r::Run) = r.records
label_of(r::Run, _="run") = r.label
kind_of(r::Run) = r.kind

# Bare vectors / pairs still work, and are assumed absolute.
records_of(x) = x isa Pair ? x.second : x
label_of(x, fallback="run") = x isa Pair ? String(x.first) : fallback
kind_of(::Any) = :absolute

"""
    greb_runs(res; ctrl="control", scnr="scenario")

Both runs from a `greb_model!` result, each tagged correctly. This is the
one-liner that keeps you out of the anomaly/absolute trap:

    runs = greb_runs(res; scnr="2xCO2")
"""
function greb_runs(res; ctrl::AbstractString="control", scnr::AbstractString="scenario")
    out = Run[]
    isempty(res.ctrl) || push!(out, control(ctrl, res.ctrl))
    isempty(res.scnr) || push!(out, scenario(scnr, res.scnr))
    out
end

"""
    monthly_climatology(records)

The final 12 months of `records`, which is the control climatology
`apply_scenario_anomalies` subtracted (see `postprocess.jl`).
"""
function monthly_climatology(records::AbstractVector)
    n = length(records)
    n < 12 && error("need at least 12 months for a climatology, got $n")
    records[(n - 11):n]
end

"""
    to_absolute(anomaly_run, control_run)

Undo `apply_scenario_anomalies`: add the control's monthly climatology back,
so a scenario can be plotted in real units alongside its control.

    abs2x = to_absolute(scenario("2xCO2", res.scnr), control("ctrl", res.ctrl))

Calendar months are matched the same way `apply_scenario_anomalies` matched
them (`mod(index-1, 12) + 1`), so this is an exact inverse.
"""
function to_absolute(anom, ctrl)
    arecs, crecs = records_of(anom), records_of(ctrl)
    kind_of(anom) === :anomaly ||
        error("to_absolute expects an :anomaly run; '$(label_of(anom))' is :$(kind_of(anom))")
    clim = monthly_climatology(crecs)
    flds = fieldnames(typeof(first(arecs)))
    recs = map(enumerate(arecs)) do (i, r)
        ref = clim[mod(i - 1, 12) + 1]
        NamedTuple{flds}(map(f -> getfield(r, f) .+ getfield(ref, f), flds))
    end
    Run(label_of(anom) * " (absolute)", recs; kind=:absolute)
end

"True if `runs` mixes absolute and anomaly quantities."
mixed_kinds(runs) = length(unique(kind_of.(runs))) > 1

"Suffix for an axis label, so a plot always says which space it is in."
kind_suffix(k::Symbol) = k === :anomaly ? " anomaly" : ""
