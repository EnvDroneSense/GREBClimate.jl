# Tutorial

This walks through the same flow as [`examples/run_greb.jl`](https://github.com/EnvDroneSense/GREBClimate.jl/blob/main/examples/run_greb.jl):
load input data, configure an experiment, run the model, and inspect the
result.

## 1. Load input data

```julia
using GREBClimate

jld2_dir = greb_data_dir()                            # see below
fields = load_greb_jld2!(jld2_dir; dataset = :ncep)   # or :era
```

[`greb_data_dir`](@ref) returns the dataset directory, downloading and caching
it (~353 MB) on first use if no local copy is found. It checks an explicit path,
then `$GREB_DATA`, then `greb_input_data/` beside the package, and only then the
network - so if you already have the data, nothing is fetched. Pass a path
directly if you prefer: `greb_data_dir("/path/to/greb_input_data")`, or hand
`load_greb_jld2!` the path itself. See [Input data](@ref).

`fields` is a [`ClimateFields`](@ref) - climatology, grid geometry, flux
corrections, and the regional-CO₂ mask/solar table. Every physics function
takes it as an explicit argument; nothing is shared as module-global state,
so you can hold several independent `fields` instances (e.g. for parameter
sweeps) in the same session.

!!! warning "`load_greb_jld2!` returns the data - it does not set globals"
    The returned `fields` must be passed to [`greb_model!`](@ref) explicitly
    (step 3). A bare [`ClimateFields`](@ref) is all zeros, and stepping the
    model on a zero climatology runs to completion and returns NaN in
    every field. `greb_model!` therefore refuses unloaded fields; see [Data-free runs](@ref) below.

## 2. Configure the experiment

[`preset`](@ref) returns the [`Config`](@ref) of a named experiment
([`preset_names`](@ref) lists them):

```julia
cfg = preset(:full_model)   # or :co2_double, :elnino, :ssp585, ...
```

A `Config` is immutable. Change a preset's physics with a NamedTuple of
[`Processes`](@ref) or [`Hydrology`](@ref) options, or build a config from its
parts; the [Configuration](@ref) page lists every option.

```julia
cfg = preset(:co2_double; hydrology = (rain = :rh,))
cfg = preset(:custom_co2; path = "my_co2.txt")                     # "year CO2" per line
cfg = preset(:decon_mean_climate; processes = (ocean = :none,))
cfg = preset(:decon_2xco2; processes = (clouds = :uniform,))
cfg = preset(:obliquity; index = 3)
cfg = preset(:earth_sun_distance; pct = 1.5)
cfg = Config(scenario = Scenario(co2 = ConstantCO2(500)))
```

## 3. Run the model

[`greb_model!`](@ref) takes a [`RunSpec`](@ref) (how many years of control
and scenario to run) and the config:

```julia
run = RunSpec(ctrl = 5, scnr = 15)
result = greb_model!(run, cfg; jld2_dir = jld2_dir, fields = fields)
```

This runs, in order: a flux-correction spin-up that holds the control climate
at the observed climatology, a control run, and a scenario run under the
experiment's forcing. The spin-up and control run at the scenario's control
CO₂, 340 ppm (280 ppm for the IPCC CO₂-table scenarios); the scenario sets its
own CO₂.

The spin-up is part of the config: `corrections = SpinUp(3)`, the original
GREB spin-up, is the default. `Stored()` uses the dataset's corrections
without a spin-up. They were computed with the MSCM physics, so with the
default hydrology the control drifts (about +2 K over 5 years in a
`:full_model` run), which shows up in the scenario anomaly. `NoCorrections()`
runs without any; with the default hydrology such a run diverges within about
15 years.

## 4. Inspect results

```julia
result.ctrl    # Vector{MonthlyRecord}, one per control-run month
result.scnr    # Vector{MonthlyRecord}, one per scenario-run month
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

A plain `mean` over the grid over-weights the polar rows; weight by
`cos(latitude)` for a global mean:

```julia
w = [cosd(-90 + (j - 0.5) * 180 / 48) for _ in 1:96, j in 1:48]
Ts_global_mean = [sum(rec.Ts .* w) / sum(w) for rec in result.ctrl]
```

The [Plots and notebook](@ref) page shows how to plot a result.

## Data-free runs

Some runs legitimately need no dataset: tests that exercise configuration or
CO₂-scenario plumbing rather than physics, and the package's own
precompilation, which must not require a 353 MB download. These opt in
explicitly:

```julia
greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = "", allow_uninitialized = true)
```

Results from such a run are structurally valid but physically meaningless -
use them to check shapes and code paths, never climate numbers.

## Next steps

- The [API Reference](@ref) lists every exported function and type.
- The [Configuration](@ref) page documents every option and preset.
- The [Model overview](@ref) explains what each component computes.
- [Plots and notebook](@ref) shows how to plot a result and explore it interactively.
- `benchmark/run_benchmarks.jl` micro-benchmarks the per-timestep physics
  kernels (`year`, `stages`, `threads`, `alloc` and `years` modes - the last
  runs a `--ctrl`/`--scnr`-year control+scenario run, for long-run cost and
  stability rather than one year).
- `test/runtests.jl` doubles as executable documentation for individual
  kernels' behavior under different configurations.
