# Tutorial

This explains the steps of a run: load input data, configure an experiment, run
the model, and inspect the result. The [Examples](examples.md) page shows the
same steps as scripts you can run, and an ensemble that asks which processes
carry the warming.

## 1. Load input data

```julia
using GREBClimate

jld2_dir = greb_data_dir()                            # see below
fields = load_climatology(jld2_dir; dataset = :ncep)   # or :era
```

[`greb_data_dir`](@ref) returns the dataset directory. It looks in this order:
an explicit path, `$GREB_DATA`, `greb_input_data/` in the package's top folder,
and an earlier download. Only when none of them exists does it download and cache
the dataset (a 321 MB download), so a local copy is never fetched again. You can pass the
path yourself, as in `greb_data_dir("/path/to/greb_input_data")`, or give
`load_climatology` the path directly. See [Input data](@ref).

`fields` is a [`ClimateFields`](@ref): the climatology, the grid geometry, the
flux corrections, and the regional-CO₂ mask and solar table. Every physics
function takes it as an argument and nothing is kept in module-global state, so
one session can hold several independent `fields`, for example for parameter
sweeps.

!!! warning "`load_climatology` returns the data and sets no globals"
    Pass the returned `fields` to [`greb_model!`](@ref) (step 3). A bare
    [`ClimateFields`](@ref) is all zeros, and a model stepped on a zero climatology
    finishes without an error and returns NaN in every field. `greb_model!`
    therefore refuses unloaded fields; see [Data-free runs](@ref) below.

## 2. Configure the experiment

[`preset`](@ref) returns the [`Config`](@ref) of a named experiment
([`preset_names`](@ref) lists them):

```julia
cfg = preset(:full_model)   # or :co2_double, :elnino, :ssp585, ...
```

A `Config` is immutable. To change a preset's physics, pass a NamedTuple of
[`Processes`](@ref) or [`Hydrology`](@ref) options, or build a config from its
parts. The [Configuration](@ref) page lists every option.

```julia
cfg = preset(:co2_double; hydrology = (rain = :rh,))
cfg = preset(:custom_co2; path = "my_co2.txt")                     # "year CO2" per line
cfg = preset(:decon_mean_climate; processes = (ocean = :none,))
cfg = preset(:decon_2xco2; processes = (clouds = :uniform,))
cfg = preset(:obliquity; index = 3)
cfg = preset(:earth_sun_distance; percent = 1.5)
cfg = Config(scenario = Scenario(co2 = ConstantCO2(500)))
```

## 3. Run the model

[`greb_model!`](@ref) takes a [`RunSpec`](@ref) (how many years of control
and scenario to run) and the config:

```julia
run = RunSpec(ctrl = 5, scnr = 15)
result = greb_model!(run, cfg; jld2_dir = jld2_dir, fields = fields)
```

A run has three phases, in this order: a flux-correction spin-up that holds the
control climate at the observed climatology, a control run, and a scenario run
under the experiment's forcing. The spin-up and the control run use the control
CO₂: 340 ppm, or 280 ppm for the presets that read a CO₂ table (`:rcp*`,
`:ssp*`, `:historical_co2`, `:a1b`, `:custom_co2`). The scenario sets its own CO₂.

The run reports its progress through the logger, one line per phase and per
simulated year:

```
[ Info: Flux-correction spin-up: CO2 = 340.0 ppm, 3 yr
[ Info: spin-up year 1: Ts global mean 13.82 °C; 178 E 9 N 26.79; 58 E 51 N 4.85
[ Info: Control run: CO2 = 340.0 ppm, 5 yr
[ Info: 1970: Ts global mean 13.82 °C; 178 E 9 N 26.79; 58 E 51 N 4.85
```

Each yearly line gives the area-weighted annual mean of the surface temperature
and its value in two cells. To silence these lines, run under a null logger:

```julia
using Logging
result = with_logger(NullLogger()) do
    greb_model!(run, cfg; jld2_dir = jld2_dir, fields = fields)
end
```

The spin-up is part of the config. The default is `corrections = SpinUp(3)`, the
original GREB spin-up. `Stored()` uses the dataset's corrections and skips the
spin-up. Those corrections were computed with the MSCM physics, so with the
default hydrology the control drifts (about +2 K over 5 years in a `:full_model`
run), and the drift shows up in the scenario anomaly. `NoCorrections()` runs
without any corrections; with the default hydrology such a run diverges within
about 15 years.

## 4. Inspect results

```julia
result.ctrl    # Vector{MonthlyRecord}, one per control-run month
result.scnr    # Vector{MonthlyRecord}, one per scenario-run month
result.ctrl_time[1]   # (year = 1970, month = 1): the time of result.ctrl[1]
result.scnr_time[1]   # (year = 1950, month = 1) for a scenario starting in 1950
result.scnr_anomaly   # true: result.scnr is the change against the control; false: absolute values
```

Each [`MonthlyRecord`](@ref) is a `NamedTuple` with fields
`Ts, Ta, To, q, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens, olr, lwdown` -
each a `(96, 48)` matrix of that month's mean. `olr` is the longwave leaving
to space and `lwdown` the longwave reaching the surface from the air.

`result.ctrl` is in absolute units. `result.scnr` is an **anomaly**: each
month minus the same calendar month of the control's final year. It stays
absolute for the orbital experiments (`:obliquity`, `:eccentricity`,
`:earth_sun_distance`, whose scenario `output` is `:absolute`) and when
`ctrl = 0`.

A plain `mean` over the grid gives the polar rows too much weight; for the
surface temperature it reads 8 to 11 K too cold. [`global_mean`](@ref) weights
each row by the cosine of its latitude:

```julia
Ts_global_mean = [global_mean(rec.Ts) for rec in result.ctrl]
```

The year and month of record `i` are `result.ctrl_time[i]` and
`result.scnr_time[i]`:

```julia
july = [rec.Ts for (rec, time) in zip(result.scnr, result.scnr_time) if time.month == 7]
```

The model's own calendar has 365 days, 730 steps a year and no leap years.
[`step_of_year`](@ref), [`day_of_year`](@ref), [`month_of_step`](@ref) and
[`decimal_year`](@ref) convert between them.

The [Plots and notebook](@ref) page shows how to plot a result.

## 5. Run many configurations

[`run_ensemble`](@ref) runs a list of configs side by side and returns one entry
per config, in the same order. Each member gives exactly the result of running
it alone. The function passed first is applied to each result inside its task,
so only what it returns is kept:

```julia
configs = [preset(:co2_double; hydrology = (evaporation = e,))
           for e in (:original, :skin, :skin_gust)]
warming = run_ensemble(RunSpec(ctrl = 1, scnr = 30), configs;
                       fields = fields, jld2_dir = jld2_dir) do result
    sum(global_mean(rec.Ts) for rec in result.scnr[end-11:end]) / 12
end
```

Start Julia with several threads (`julia -t 8,0`) so that the members run in
parallel. `ntasks` sets how many run at once (default: the number of threads).
Each task works on its own copy of `fields`, about 370 MB. A member whose
physics differs from the default needs its own spin-up (the default
`SpinUp(3)`), because the stored corrections fit one configuration only.

## 6. Check a run

`greb_model!` takes an `observer`, a function that is called twice per step with
the state. Two observers are included. `GREBClimate.RangeCheck()` catches a run
that leaves the physical range without turning non-finite, which could otherwise
pass for a result:

```julia
check = GREBClimate.RangeCheck()
result = greb_model!(run, cfg; jld2_dir = jld2_dir, fields = fields, observer = check)
GREBClimate.in_range(check)   # false if Ts, Ta, To or q left its range
check.first                   # the first step outside, or nothing
check.seen                    # lowest and highest value of each field
```

`GREBClimate.BudgetCheck()` checks the model's bookkeeping: each store must
change by the sum of its flows. Both have docstrings in the
[API Reference](@ref), and [Diagnostics](dev/diagnostics.md) shows a budget report.

## Data-free runs

Some runs need no dataset: tests of the configuration or the CO₂-scenario
plumbing rather than the physics, and the package's own precompilation, which
must not require a 321 MB download. These runs opt in explicitly:

```julia
greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = "", allow_uninitialized = true)
```

The results of such a run have the right shape but no physical meaning. Use them
to check shapes and code paths, never to read climate numbers.

## Next steps

- The [API Reference](@ref) lists every exported function and type.
- The [Configuration](@ref) page documents every option and preset.
- The [Model overview](@ref) explains what each component computes.
- [Plots and notebook](@ref) shows how to plot a result and explore it interactively.
- The `test/` directory has one file per subject and shows how the kernels behave under different configurations.
