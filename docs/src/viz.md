# Plots and notebook

The repository ships a plotting toolbox (`viz/`) and a Pluto notebook
(`notebooks/GREB_explorer.jl`) built on it. Neither is part of the package:
they share their own environment, `viz/Project.toml` (Plots, Pluto, PlutoUI).
Both need the dataset on disk (see [Input data](@ref)); they never download it.

## The notebook

From the repository root:

```bash
julia notebooks/launch_pluto.jl
```

This opens Pluto in the browser with the explorer loaded; the first launch also
sets up the `viz/` environment, which takes a few minutes. Pick an experiment
and the number of control and scenario years, tick **Run the model**, then
choose a variable. The notebook runs a 3-year flux-correction spin-up before
the control, and shows the plots below plus a frame-by-frame view of the run
that can be played or exported as a GIF.

## From a script

Set up the environment once (the notebook launcher does this for you):

```bash
julia viz/setup.jl
```

Then, with `julia --project=viz`:

```julia
include("viz/GREBViz.jl"); using .GREBViz

plot_map(result; var = :Ts, month = :mean, fields = fields)   # fields adds coastlines
plot_timeseries(result; annual = true)
plot_seasonal(result; var = :precip)
plot_hovmoller(result)

ev = evolution(result; var = :Ts, step = :year)
evolution_frame(ev, 3; fields = fields)
evolution_gif("run.gif", ev; fields = fields)
```

Every plot takes a `greb_model!` result, drawing control and scenario as
separate panels, or a plain vector of monthly records. The scenario is usually
an anomaly (a change against the control), drawn on a colour scale centred at
zero: red is warmer or more, blue colder or less.

## What each plot shows

| Function | Shows | How to read it |
|:---------|:------|:---------------|
| `plot_map` | One variable over the globe, averaged over the run (`month = :mean`), in the last month (`:last`) or in month `n` | *Where* a quantity is high or low, or, for the scenario, where it changed most |
| `plot_timeseries` | The global mean per month, or per year with `annual = true` | *How fast* the climate changes, and whether it has settled |
| `plot_seasonal` | The average of each calendar month over the whole years in the run | The seasonal cycle, and for the scenario which seasons changed most |
| `plot_hovmoller` | The mean over longitude, as a colour against month (x) and latitude (y) | How a signal moves between latitudes and seasons, e.g. stronger polar warming each autumn |
| `evolution` + `evolution_frame` / `evolution_gif` | The map month by month (or year by year) on a fixed colour scale | The map changing through the run; same colour means same value in every frame |

The numbers behind the plots are available directly: `series`, `annual`,
`seasonal_cycle`, `field`, `hovmoller`, `map_frames` and `area_weights`.

```bash
julia --project=viz viz/demo.jl [output_dir]   # runs GREB and saves every plot
julia --project=viz viz/test.jl                # the toolbox's tests
```
