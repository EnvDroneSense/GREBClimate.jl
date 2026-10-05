# GREB Climate Model - Julia Implementation

<img src="docs/src/assets/logo-wordmark.svg" width="400" alt="GREBClimate.jl logo">

[![CI](https://github.com/EnvDroneSense/GREBClimate.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/EnvDroneSense/GREBClimate.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/github/EnvDroneSense/GREBClimate.jl/graph/badge.svg?token=CKFBW810SH)](https://codecov.io/github/EnvDroneSense/GREBClimate.jl)
[![docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://EnvDroneSense.github.io/GREBClimate.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.10+-9558B2?logo=julia)](https://julialang.org/)
[![Pluto](https://img.shields.io/badge/Pluto-Interactive-purple)](https://github.com/fonsp/Pluto.jl)

A high-performance Julia translation of the **Globally Resolved Energy Balance (GREB)** climate model, originally developed by Dietmar Dommenget and colleagues at Monash University.

GREBClimate is a **Julia package**: you call it from a script or the REPL, and nothing is held as module-global state. A plotting toolbox (`viz/`) and an interactive [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebook ship alongside it. The full guide is in the **[documentation](https://EnvDroneSense.github.io/GREBClimate.jl/dev/)**.

## Table of Contents

- [About the Model](#about-the-model)
- [Features](#features)
- [Quick Start](#-quick-start)
- [Input Data](#-input-data)
- [Running the Model](#-running-the-model)
- [Model Components](#-model-components)
- [Contributing](#contributing)
- [References](#-references)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## About the Model

GREB is a conceptual climate model that simulates the global energy balance on a **3.75° × 3.75°** grid (96 longitudes × 48 latitudes). It uses a **12-hour time step** with **30-minute sub-steps** for atmospheric transport (730 steps per simulated year). The atmospheric circulation, clouds and ocean mixed-layer depth are taken from observed climatology, which keeps it fast enough for long runs and many experiments. This package is a translation of the original Fortran90 code.

## Features

- **40+ experiments**: CO₂ scaling, IPCC RCP and SSP scenarios, a historical CO₂ hindcast, your own CO₂ trajectory, solar, orbital and paleoclimate forcing, and ENSO and regional-CO₂ runs ([experiment list](https://EnvDroneSense.github.io/GREBClimate.jl/dev/configuration/#Presets))
- **Deconstruction experiments**: switch individual feedback processes off to isolate their role in the mean climate or the 2×CO₂ response ([switches](https://EnvDroneSense.github.io/GREBClimate.jl/dev/configuration/))
- **Two climatologies**: NCEP and ERA-Interim
- **Ensembles**: run many configurations side by side with `run_ensemble`, each member identical to its own single run
- **Fast**: SIMD-vectorised physics, and the temperature and humidity transport run concurrently with `julia -t 2,0`
- **Plots and notebook**: maps, time series, seasonal cycles, Hovmöller diagrams and animations ([guide](https://EnvDroneSense.github.io/GREBClimate.jl/dev/viz/))

## 🚀 Quick Start

### Prerequisites

Requires **Julia 1.10** (the current LTS) or later. Download from [julialang.org](https://julialang.org/downloads/).

### Installation

GREBClimate is registered in the Julia General Registry:

```julia
using Pkg
Pkg.add("GREBClimate")
```

To work on the package itself, or to use the notebook below, clone the repository and instantiate it:

```bash
git clone https://github.com/EnvDroneSense/GREBClimate.jl
cd GREBClimate.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Upgrading from 1.x

Version 2.0 replaces `PhysicsConfig` and `create_experiment_config` with `preset(name)` and `Config`, moves the spin-up length from `RunSpec` into the config (`SpinUp(years)`), and renames some presets. A few experiments give different results. [CHANGELOG.md](CHANGELOG.md) lists every change.

### Launch the notebook (optional)

```bash
julia notebooks/launch_pluto.jl
```

The first launch sets up its `viz/` environment. See [Plots and notebook](https://EnvDroneSense.github.io/GREBClimate.jl/dev/viz/) for the plots and scripting.

## 📂 Input Data

The model reads a ~439 MB JLD2 dataset of climatologies, flux corrections and forcing tables. It is not in the repository: `greb_data_dir()` finds a local copy (an explicit path, `$GREB_DATA`, or `greb_input_data/` beside the repository) and otherwise downloads and caches it once (~353 MB). Set `DATADEPS_ALWAYS_ACCEPT=true` in non-interactive sessions. The layout, the loading rules and how to regenerate the data are in [Input data](https://EnvDroneSense.github.io/GREBClimate.jl/dev/data/).

## 🎮 Running the Model

```julia
using GREBClimate

dir    = greb_data_dir()
fields = load_climatology(dir; dataset=:ncep)            # returns the data; pass it on
cfg    = preset(:co2_double)                            # preset_names() lists them all
result = greb_model!(RunSpec(ctrl=5, scnr=15), cfg; jld2_dir=dir, fields=fields)
```

A run has three phases, in years:

| Phase | What it does |
|:------|:-------------|
| Spin-up | Derives the flux corrections holding the control at the observed climate: `corrections = SpinUp(3)` in the config, the default. `Stored()` uses the dataset's corrections, `NoCorrections()` lets the control drift. |
| `ctrl` (`RunSpec`) | Control run at 340 ppm CO₂ (280 ppm for the IPCC CO₂-table scenarios; `:rcp85_boundary` uses 340) |
| `scnr` (`RunSpec`) | Scenario run under the experiment's forcing |

`result.ctrl` and `result.scnr` are vectors of monthly means, each a `NamedTuple` of 96×48 fields (`Ts, Ta, To, q, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens, olr, lwdown`). `result.scnr` is an **anomaly** against the control's final year, except for the orbital experiments and runs with `ctrl=0`. `global_mean(rec.Ts)` gives an area-weighted global mean, and `run_ensemble` runs a list of configurations side by side. The [Tutorial](https://EnvDroneSense.github.io/GREBClimate.jl/dev/tutorial/) covers configuration switches, experiment keywords, reading the results, ensembles and run checks; [`examples/run_greb.jl`](examples/run_greb.jl) is a runnable script.

## 🔬 Model Components

| Component | What it does |
|:----------|:-------------|
| Shortwave radiation | Absorbed sunlight, with ice-albedo feedback and climatological clouds |
| Longwave radiation | Emission and back-radiation; emissivity from CO₂, water vapor and clouds (the greenhouse effect) |
| Hydrology | Evaporation, precipitation and latent heat |
| Atmospheric transport | Diffusion and advection of heat and moisture by climatological winds (~93% of run time) |
| Ocean | Mixed-layer heat content, exchange with the deep ocean, and sea ice |
| Flux corrections | Keep the control climate at the observed climatology |

The [Model overview](https://EnvDroneSense.github.io/GREBClimate.jl/dev/model/) gives the energy-balance equation, the functions behind each component, and how a run is structured.

## Contributing

Bug reports and contributions are welcome. [CONTRIBUTING.md](CONTRIBUTING.md) covers the development setup, the tests, the project layout, how to report a bug and the roadmap.

## 📚 References

### Primary Publications

1. **Dommenget, D., and Flöter, J. (2011)**. Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model. *Journal of Climate Dynamics*, 37: 2143. [doi:10.1007/s00382-011-1026-0](https://doi.org/10.1007/s00382-011-1026-0)

2. **Stassen, C., Dommenget, D., and Loveday, N. (2019)**. A hydrological cycle model for the Globally Resolved Energy Balance (GREB) model v1.0. *Geoscientific Model Development*, 12, 425-440. [doi:10.5194/gmd-12-425-2019](https://doi.org/10.5194/gmd-12-425-2019)

3. **Dommenget, D., Nice, K., Bayr, T., Kasang, D., Stassen, C., and Rezny, M.** The Monash Simple Climate Model Experiments: An interactive database of the mean climate, climate change and scenarios simulations. *Geoscientific Model Development*, 12, 2155-2179. [doi:10.5194/gmd-12-2155-2019](https://doi.org/10.5194/gmd-12-2155-2019)

### Original GREB Model
- [Monash University GREB Homepage](https://users.monash.edu.au/~dietmard/content/GREB/code.html)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Original GREB Model**: Dietmar Dommenget, Janine Flöter, Tobias Bayr, Christian Stassen (Monash University)
- **Julia Translation and Optimization**: Thomas Struys (UGent)
- **Julia Development Guidance and Initial Package Refactor**: Michiel Stock (UGent)
- **Pluto.jl**: For the interactive notebook environment
- **Julia Community**: For excellent scientific computing tools

---

**Contact**: For questions or support, please open an issue on GitHub.
