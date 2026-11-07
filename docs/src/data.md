# Input data

The model reads a **JLD2** dataset ([JuliaIO/JLD2.jl](https://github.com/JuliaIO/JLD2.jl)):
climatology, flux corrections, solar forcing and scenario tables, under one
directory. Each field file stores an `Array{Float32}` under `"data"`, plus
`"dim_names"`, `"description"`, `"units"`, `"source"`, `"period"` and
`"dataset_version"`, and optionally `"coords"` (longitude and latitude in
degrees, time in days from 1 January) and `"ctl"` (the original GrADS
metadata). [`read_field`](@ref) returns the extra entries as `info`. The model computes in `Float32` too, so nothing is converted on load.

## Getting the data

At ~439 MB unpacked (a 321 MB download) the dataset is not committed to the repository.
[`greb_data_dir`](@ref) fetches and caches it on first use via
[DataDeps.jl](https://github.com/oxinabox/DataDeps.jl):

```julia
using GREBClimate
dir    = greb_data_dir()        # asks, downloads (321 MB) and caches on first use
fields = load_climatology(dir; dataset = :ncep)   # or :era
```

It resolves in this order; only the last step touches the network:

| | Source |
|--:|:-------|
| 1 | an explicit path, `greb_data_dir("/path/to/greb_input_data")` |
| 2 | `ENV["GREB_DATA"]` |
| 3 | `greb_input_data/` in the package's top folder |
| 4 | the `GREB-input-data-v2` DataDep, if it was downloaded before |
| 5 | the same DataDep, downloaded (after asking) and SHA256-verified |

The download is cached under `~/.julia/scratchspaces/.../datadeps/GREB-input-data-v2`,
so it happens once per machine. A dataset downloaded for GREBClimate 2.2.x
(`GREB-input-data`) is not used and not deleted.

### Packs

The dataset comes in five packs. The first download asks once whether to fetch
the whole dataset (321 MB). Yes fetches all five packs without further questions.
No fetches the core pack only, which is enough for a default run; `packs = :core`
skips the question. A run that needs an optional pack stops with an error that
names it, and [`download_greb_pack`](@ref) asks for and adds it. A cached dataset
or a folder of your own is never extended automatically.

| Pack | Files | Download | Needed for |
|:-----|------:|---------:|:-----------|
| `core` | 17 | 103 MB | every run |
| `era` | 4 | 38 MB | `load_climatology(dir; dataset = :era)` |
| `cmip5` | 5 | 54 MB | the `:rcp85` experiment |
| `enso` | 10 | 121 MB | the `:elnino` and `:lanina` experiments |
| `solar` | 3 | 4 MB | the `:paleo_231kyr`, `:eccentricity` and `:obliquity` experiments |

Two environment variables matter:

| Variable | Effect |
|:---------|:-------|
| `DATADEPS_ALWAYS_ACCEPT=true` | Skip the download prompts and fetch everything. **Required** in CI or any non-interactive session, which otherwise blocks waiting on stdin. |
| `DATADEPS_DISABLE_DOWNLOAD=true` | Make a would-be download throw instead. |

The tests and benchmarks never download: they resolve with
`allow_download=false`, which stops after step 4, and skip data-dependent work
when no local copy exists.

In Pluto there is no terminal to answer the question. Run
`using GREBClimate; greb_data_dir()` once in a REPL first; the download is cached
for every project.

## Loading

[`load_climatology`](@ref) *returns* a [`ClimateFields`](@ref) holding the
climatology, grid geometry, flux corrections and solar table. It does not set
global state: pass the value to [`greb_model!`](@ref) with `fields = ...`. You
can hold several independent instances in one session, e.g. for parameter
sweeps.

A bare `ClimateFields()` is all zeros, so `greb_model!` refuses it; runs that need no data opt in explicitly (see [Data-free runs](@ref)).

## Directory structure

```
greb_input_data/                    # 39 files, ~439 MB
├── static/
│   ├── global.topography.jld2      # 2D (96×48)
│   └── greb.glaciers.jld2          # 2D (96×48)
├── climatology/                    # 31 files; all 3D (96×48×730) unless noted
│   │   # dataset=:ncep
│   ├── ncep.tsurf.1948-2007.clim.jld2
│   ├── ncep.zonal_wind.850hpa.clim.jld2
│   ├── ncep.meridional_wind.850hpa.clim.jld2
│   ├── ncep.atmospheric_humidity.clim.jld2
│   ├── ncep.soil_moisture.clim.jld2        # also used by dataset=:era
│   │   # dataset=:era (alternative to the ncep.* four above)
│   ├── erainterim.tsurf.1979-2015.clim.jld2
│   ├── erainterim.zonal_wind.850hpa.clim.jld2
│   ├── erainterim.meridional_wind.850hpa.clim.jld2
│   ├── erainterim.atmospheric_humidity.clim.jld2
│   │   # common to both datasets
│   ├── isccp.cloud_cover.clim.jld2
│   ├── woce.ocean_mixed_layer_depth.clim.jld2
│   ├── Tocean.clim.jld2
│   ├── erainterim.omega.vertmean.clim.jld2
│   ├── erainterim.omega_std.vertmean.clim.jld2
│   ├── erainterim.windspeed.850hpa.clim.jld2
│   ├── flux_corrections.jld2       # Tsurf/vapor/Tocean corrections
│   │   # CMIP5 RCP8.5 anomalies - climate-change experiments only
│   ├── cmip5.{tsurf,zonal.wind,meridional.wind,windspeed,omega}.rcp85.ensmean.forcing.jld2
│   │   # ENSO anomalies - :elnino / :lanina only (10 files)
│   └── erainterim.{tsurf,zonal.wind,meridional.wind,windspeed,omega}.{elnino,lanina}.forcing.jld2
├── solar/
│   └── solar_radiation.clim.jld2   # 2D (48×730)
├── solar_scenarios/                # pack solar
│   ├── solar_paleo.jld2
│   ├── solar_eccentricity.jld2     # (ecc_index, lat, time)
│   └── solar_obliquity.jld2        # (obl_index, lat, time)
└── scenario/                       # core
    ├── ipcc_scenarios.jld2         # "scenarios", "labels", "periods"; keyed "rcp26" ... "rcp85", "ssp119" ... "ssp585", "hist"
    └── historical_emissions_population.jld2   # not read by the model (11 KB)
```

The CO2 tables of `ipcc_scenarios.jld2` carry a label, read with
[`co2_scenario_label`](@ref): the `rcp*` tables hold an equivalent concentration
(`:co2_equivalent`), the `ssp*` tables and `hist` hold CO2 only (`:co2`). The SSP
tables use the historical record before 2009.

## Regenerating the data (maintainers)

Not needed to run the model. The raw GREB `.bin` inputs are collated from
several upstream sources and are not redistributed;
[`DATA_README.md`](https://github.com/EnvDroneSense/GREBClimate.jl/blob/main/DATA_README.md)
documents them.

```bash
julia --project=. tools/dataset/convert_greb_to_jld2.jl <input_dir> [output_dir]   # default: greb_input_data/
julia --project=. tools/dataset/package_dataset.jl greb_input_data dataset_packs
```

The second command checks the tree against the tools' file list, builds the five
reproducible `tar.xz` archives (needs GNU `tar` and `xz`) and prints the
`DATA_PACKS` block to paste into `src/data.jl` before attaching the archives to
the `data-v2` release. The full reference is the developer page "Dataset tools".
