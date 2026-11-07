# Plots and notebook

The repository ships a plotting toolbox (`viz/`, module `GREBViz`) and a Pluto notebook
(`viz/notebooks/GREB_explorer.jl`) built on it. Neither is part of the package:
they share their own environment, `viz/Project.toml` (Plots, Pluto, PlutoUI).
Both need the dataset on disk (see [Input data](@ref)); they never download it.

| Page | Holds |
|:-----|:------|
| This page | Setting up, the notebook, running the toolbox from a script |
| [The plots](plots.md) | Every plot with an example call and what it shows |
| [Building blocks](tools.md) | The functions under the plots: reductions, masks, outlines. For making your own plots |
| [Reading the plots](interpreting.md) | Colour scales, anomalies, area weighting and the limits of a 96x48 grid |

## The notebook

From the repository root:

```bash
julia viz/notebooks/launch_pluto.jl
```

This opens Pluto in the browser with the explorer loaded. Each launch brings
the `viz/` environment up to date; the first takes a few minutes. Pick an experiment
and the number of control and scenario years, tick **Run the model**, then
choose a variable. The notebook runs the preset's flux-correction spin-up
(3 years for most presets) before the control, shows the plots of the next page, and adds a
frame-by-frame view of the run that can be played or exported as a GIF.

## From a script

Set up the environment (the notebook launcher does this on every start). Run it
again after updating the repository:

```bash
julia viz/setup.jl
```

Then, with `julia --project=viz`:

```julia
include("viz/GREBViz.jl"); using .GREBViz

plot_map(result; var = :Ts, fields = fields)   # fields adds coastlines
plot_timeseries(result; annual = true)
```

`result` is what `greb_model!` returns and `fields` the climatology from
`load_climatology`. Two commands run the whole toolbox:

```bash
julia --project=viz viz/demo.jl [output_dir]   # runs GREB and saves every plot
julia --project=viz viz/test.jl                # the toolbox's tests, no dataset needed
julia --project=viz -t 2,0 viz/docs_figures.jl # regenerates the figures of these pages
```
