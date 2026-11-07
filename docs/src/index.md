# GREBClimate.jl

A Julia translation of the **Globally Resolved Energy Balance (GREB)**
climate model, originally developed by Dietmar Dommenget and colleagues at
Monash University.

It runs on a fixed 96 × 48 grid with a 12-hour step; the
[Model overview](@ref) explains its components.

## Installation

Requires Julia 1.10 or later.

```julia
using Pkg
Pkg.add("GREBClimate")
```

Or, working from a clone:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The model needs a ~439 MB input dataset, which [`greb_data_dir`](@ref)
downloads and caches on first use; see [Input data](@ref).

## Where to go next

| Page | For |
|:-----|:----|
| [Tutorial](@ref) | Load the data, configure an experiment, run it and read the result |
| [Input data](@ref) | Where the dataset comes from, how it is found, and its layout |
| [Model overview](@ref) | What the model computes and how a run is structured |
| [Time and calendar](@ref) | The model calendar, the time value of a run, and what the time system does not do yet |
| [Plots and notebook](@ref) | The `viz/` plotting toolbox and the Pluto explorer |
| [Configuration](@ref) | The parts of a configuration, their options, and the experiment presets |
| [API Reference](@ref) | Every exported function and type |
