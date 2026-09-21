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
    model on a zero climatology runs to completion while producing a
    physically meaningless world pinned at the 40 K stability floor
    (−233 °C). `greb_model!` therefore
    refuses unloaded fields; see [Data-free runs](@ref) below.

## 2. Configure the experiment

[`create_experiment_config`](@ref) returns a [`PhysicsConfig`](@ref) preset
for a named experiment:

```julia
cfg = create_experiment_config(:full_model)   # or :co2_double, :elnino, :rcp85, ...
```

`cfg` is a mutable struct - override individual switches after construction,
e.g. `cfg.log_rain = 1` to pick a different hydrology parameterization.
See the [Physics Switches](@ref) page for the full list of switches and
what each one controls.

Some experiments take extra keywords:

```julia
cfg = create_experiment_config(:custom_co2; co2_path = "my_co2.txt")    # "year CO2" per line
cfg = create_experiment_config(:decon_mean_climate; log_ocean_dmc = false)
cfg = create_experiment_config(:decon_2xco2; log_clouds_drsp = false)
cfg = create_experiment_config(:obliquity; orbital_index = 3)
cfg = create_experiment_config(:earth_sun_distance; earth_sun_distance_pct = 1.5)
```

The `log_*` keywords apply only to the two `:decon_*` experiments; passing one
elsewhere warns.

## 3. Run the model

[`greb_model!`](@ref) takes a [`RunSpec`](@ref) (how many years of
flux-correction spin-up, control, and scenario to run) and the config:

```julia
run = RunSpec(flux = 3, ctrl = 5, scnr = 15)
result = greb_model!(run, cfg; jld2_dir = jld2_dir, fields = fields)
```

This runs, in order: a flux-correction spin-up that holds the control climate
at the observed climatology, a control run, and a scenario run under the
experiment's forcing. The spin-up and control run at 340 ppm CO₂ (280 ppm for
the IPCC-scenario experiments); the experiment sets only the scenario's CO₂.

Without the spin-up (`flux = 0`) the control drifts: about +2 K over 5 years
in a `:full_model` run, which then shows up in the scenario anomaly. Use
`flux = 3`, the original GREB default, for experiments.

## 4. Inspect results

```julia
result.ctrl    # Vector{MonthlyRecord}, one per control-run month
result.scnr    # Vector{MonthlyRecord}, one per scenario-run month
```

Each [`MonthlyRecord`](@ref) is a `NamedTuple` with fields
`Ts, Ta, To, q, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens` - each
a `(96, 48)` matrix of that month's mean.

`result.ctrl` is in absolute units. `result.scnr` is an **anomaly**: each
month minus the same calendar month of the control's final year. It stays
absolute for the orbital experiments (`:obliquity`, `:eccentricity`,
`:earth_sun_distance`) and when `ctrl = 0`.

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
- The [Physics Switches](@ref) page documents every `PhysicsConfig` field.
- The [Model overview](@ref) explains what each component computes.
- [Plots and notebook](@ref) shows how to plot a result and explore it interactively.
- `benchmark/run_benchmarks.jl` micro-benchmarks the per-timestep physics
  kernels (`year`, `stages`, `threads` and `alloc` modes).
- `test/runtests.jl` doubles as executable documentation for individual
  kernels' behavior under different config switches.
