# Shared fixtures and helpers for the test suite.
#
# `allow_download=false` so running the tests can never trigger the 353 MB
# DataDeps download; data-dependent testsets @test_skip when this is absent.
const DATA_DIR = something(greb_data_dir(; allow_download=false),
    joinpath(@__DIR__, "..", "greb_input_data"))

const X, Y, N = xdim, ydim, nstep_yr

"""Run `f` with the logger muted - the model logs a progress line per year."""
quiet(f) = with_logger(f, NullLogger())

"A `Float32` grid map with the value `v` in every cell."
field(v) = fill(Float32(v), X, Y)

"""Run `f(dir)` in a fresh temp directory, removed afterwards even on failure."""
function with_tempdir(f)
    dir = mktempdir()
    try
        f(dir)
    finally
        rm(dir; recursive=true, force=true)
    end
end

"""
Run `f()` when the dataset is present; otherwise record a skipped test. CI has
no dataset, so what `f` checks runs on a developer machine only.
"""
with_dataset(f) = isdir(DATA_DIR) ? f() : @test_skip "greb_input_data/ not present"

"""
Write a field file at `path`, creating its directory: `fill(value, dims...)`
with the names of its axes, by default those of a climatology.
"""
function write_field(path, value, dims=(X, Y, N); dim_names=["lon", "lat", "time"][1:length(dims)])
    mkpath(dirname(path))
    jldopen(path, "w") do file
        file["data"] = fill(value, dims...)
        file["dim_names"] = dim_names
        file["dataset_version"] = GREBClimate.DATASET_VERSION
    end
    return path
end

"""
Write `climatology/flux_corrections.jld2` into `dir`, each table filled with
its value. A table given as `nothing` is left out of the file.
"""
function write_flux_corrections(dir; Ts=0.0f0, q=0.0f0, To=0.0f0)
    mkpath(joinpath(dir, "climatology"))
    jldopen(joinpath(dir, "climatology", "flux_corrections.jld2"), "w") do file
        for (key, value) in (("Tsurf_flux_correction", Ts), ("vapour_flux_correction", q), ("Tocean_flux_correction", To))
            value === nothing || (file[key] = fill(value, X, Y, N))
        end
        file["dataset_version"] = GREBClimate.DATASET_VERSION
    end
    return dir
end

"""
Write a minimal `scenario/ipcc_scenarios.jld2` into `dir`. `table` maps the
scenario key to a year=>ppm Dict; the label of each table is `co2_equivalent`
for an `rcp` key and `co2` for any other.
"""
function write_ipcc_scenarios(dir, table)
    mkpath(joinpath(dir, "scenario"))
    jldopen(joinpath(dir, "scenario", "ipcc_scenarios.jld2"), "w") do file
        file["scenarios"] = table
        file["labels"] = Dict(key => startswith(key, "rcp") ? "co2_equivalent" : "co2" for key in keys(table))
        file["dataset_version"] = GREBClimate.DATASET_VERSION
    end
    return dir
end

"""
Write the three `solar_scenarios/*.jld2` tables the paleo/orbital experiments
swap in at scenario start, filled with a recognisable sentinel.
"""
function write_solar_scenarios(dir; sentinel=999.0)
    write_field(joinpath(dir, "solar_scenarios", "solar_paleo.jld2"), sentinel, (Y, N); dim_names=["lat", "time"])
    for which in ("obliquity", "eccentricity")
        jldopen(joinpath(dir, "solar_scenarios", "solar_" * which * ".jld2"), "w") do file
            file["data"] = fill(sentinel, 1, Y, N)
            file["dim_names"] = ["index", "lat", "time"]
            file["coords"] = Dict(1 => [0.0])
        end
    end
    return dir
end

"""
Write the boundary-anomaly files of the three sources into `dir/climatology`,
each array filled with its own value, and return those values by source and
array. The file names are written out here, not taken from the package.
"""
function write_boundary_anomalies(dir)
    cmip5 = (Ts_clim=2.0, u_clim=3.0, v_clim=4.0, omega_clim=5.0, wind_speed_clim=6.0)
    era = (Ts_clim=7.0, u_clim=8.0, v_clim=9.0, omega_clim=10.0, wind_speed_clim=11.0)
    anomalies = Dict(:cmip5_rcp85 => cmip5, :elnino => era, :lanina => era)
    part = (Ts_clim="tsurf", u_clim="zonal.wind", v_clim="meridional.wind", omega_clim="omega",
        wind_speed_clim="windspeed")
    file(source, name) = source === :cmip5_rcp85 ? "cmip5.$(part[name]).rcp85.ensmean.forcing.jld2" :
                         "erainterim.$(part[name]).$source.forcing.jld2"
    for (source, v) in anomalies, name in keys(part)

        write_field(joinpath(dir, "climatology", file(source, name)), v[name])
    end
    return anomalies
end

"""
A `ClimateFields` non-degenerate enough to exercise the full physics pipeline
without the real dataset: finite mixed-layer depth, land in the western half,
and non-zero winds/humidity so circulation actually transports something.
"""
function synthetic_fields()
    f = ClimateFields()
    f.mld_clim .= 50.0f0
    f.z_topo[1:(X-48), :] .= 100.0f0
    for k in 1:N, j in 1:Y, i in 1:X
        f.Ts_clim[i, j, k] = 288.0f0 - 40.0f0 * abs(j - Y / 2) / (Y / 2)
        f.u_clim[i, j, k] = 5.0f0 * sinpi(2 * j / Y)
        f.v_clim[i, j, k] = 2.0f0 * cospi(2 * i / X)
        f.q_clim[i, j, k] = 0.005f0
        f.wind_speed_clim[i, j, k] = 6.0f0
    end
    f.To_clim .= 283.0f0
    f.cloud_clim .= 0.5f0
    f.soil_wetness_clim .= 0.4f0
    GREBClimate.split_winds!(f)
    return f
end

gmean(x) = sum(x) / length(x)

"""
A `ClimateFields` that is the same in every cell and at every step: topography
`z_topo` everywhere (above 0 m is land), a 50 m mixed layer, and the given
soil wetness, winds and vertical velocity. For checking a kernel against a
hand calculation.
"""
function constant_fields(; z_topo, swet=1.0, u=0.0, v=0.0, omega=0.0, omega_std=0.0, wind_speed=0.0)
    f = ClimateFields()
    f.z_topo .= z_topo
    f.mld_clim .= 50.0
    f.Ts_clim .= 280.0
    f.To_clim .= 285.0
    f.q_clim .= 0.006
    f.cloud_clim .= 0.5
    f.soil_wetness_clim .= swet
    f.u_clim .= u
    f.v_clim .= v
    f.omega_clim .= omega
    f.omega_std_clim .= omega_std
    f.wind_speed_clim .= wind_speed
    return f
end

"A `MonthlyRecord` with every field filled with `v`; a keyword sets one field to another value."
uniform_record(v; kw...) = MonthlyRecord(map(n -> fill(Float32(get(kw, n, v)), X, Y), fieldnames(MonthlyRecord)))

"The time of scenario step `step` of a scenario that started in `start`."
scenario_time(start, step=1) = ModelTime(scenario, start, step)

"The time of control step `step`."
control_time(step=1) = ModelTime(control, 1970, step)

"The times of `n` monthly records from January of `start`: record `i` is month `mod1(i, 12)` of year `start + (i - 1) ÷ 12`."
record_times(start, n) = [(year=start + (i - 1) ÷ 12, month=mod1(i, 12)) for i in 1:n]

"A fresh copy of the four state arrays `init_model!` returned in `ini`."
initial_state(ini) = (Ts=copy(ini.surf.Ts), Ta=copy(ini.surf.Ta), q=copy(ini.surf.q), To=copy(ini.surf.To))

"""
Run the model step at time `t` on the state `s` (as from `initial_state`), muted.
`state`, `ws` and `acc` are fresh unless given; pass them to carry them over to
the next step. Returns the records and record times the step wrote.
"""
function step!(t, CO2, s, fields, r; state=ModelState(), ws=ModelWorkspace(),
        acc=MonthlyAccumulator(), observer=nothing)
    records, times = MonthlyRecord[], RecordTime[]
    quiet() do
        time_loop!(t, CO2, s.Ts, s.Ta, s.To, s.q, records, times, fields, state, ws, acc, r; observer)
    end
    return (; records, times)
end

"""
Run `greb_model!(run, config; kwargs...)` with the logger muted, on `fields`
(default `synthetic_fields()`) and without the dataset: `jld2_dir` is empty
unless given, and `allow_uninitialized` is set.
"""
function run_synthetic(run, config; fields=synthetic_fields(), jld2_dir="", kwargs...)
    return quiet() do
        greb_model!(run, config; jld2_dir, fields, allow_uninitialized=true, kwargs...)
    end
end

struct StopRun <: Exception end

"""
Run `greb_model!(run, config; kwargs...)` up to the first step of `phase`
(`control` or `scenario`), return `f(view)` of the observer's view there, and stop
the run. `point` is `:after_tendencies` (state before the step's update) or
`:after_step`. For what is already decided at the first step - the CO2, the
solar table, the climatology in use - without paying for the rest of the year.
The observer is not called during the spin-up, so pass `NoCorrections()` or
`SpinUp(0)` unless the spin-up is what is being tested.
"""
function at_first_step(f, run, config; phase=scenario, point=:after_step, kwargs...)
    seen = nothing
    function observer(pt, view)
        (view.phase === phase && pt === point) || return nothing
        seen = f(view)
        throw(StopRun())
    end
    try
        quiet() do
            greb_model!(run, config; allow_uninitialized=true, observer, kwargs...)
        end
    catch e
        e isa StopRun || rethrow()
    end
    return seen
end
