# GREBClimate.jl

A Julia translation of the **Globally Resolved Energy Balance (GREB)**
climate model, originally developed by Dietmar Dommenget and colleagues at
Monash University.

The model simulates the global energy balance on a **3.75° × 3.75°** grid
(96 longitudes × 48 latitudes), stepping every 12 hours with 30-minute
sub-steps for atmospheric circulation (730 time steps per year). The
[Model overview](@ref) explains its components.

## Installation

Requires Julia 1.10 (the current LTS) or later.

```julia
using Pkg
Pkg.develop(url = "https://github.com/EnvDroneSense/GREBClimate.jl")
```

Or, working from a clone:

```julia
julia --project=.
using Pkg; Pkg.instantiate()
```

The model needs a ~439 MB input dataset, which [`greb_data_dir`](@ref)
downloads and caches on first use; see [Input data](@ref).

## Where to go next

| Page | For |
|:-----|:----|
| [Tutorial](@ref) | Load the data, configure an experiment, run it and read the result |
| [Input data](@ref) | Where the dataset comes from, how it is found, and its layout |
| [Model overview](@ref) | What the model computes and how a run is structured |
| [Plots and notebook](@ref) | The `viz/` plotting toolbox and the Pluto explorer |
| [Physics Switches](@ref) | Every `PhysicsConfig` field and experiment preset |
| [API Reference](@ref) | Every exported function and type |
