# Input data

The model reads a **JLD2** dataset ([JuliaIO/JLD2.jl](https://github.com/JuliaIO/JLD2.jl)):
climatology, flux corrections, solar forcing and scenario tables, under one
directory. Each field file stores an `Array{Float32}` under `"data"`, plus
`"dim_names"` and optionally `"coords"` and `"ctl"` (the original GrADS
metadata). The model computes in `Float32` too, so nothing is converted on load.

## Getting the data

At ~439 MB unpacked the dataset is not committed to the repository.
[`greb_data_dir`](@ref) fetches and caches it on first use via
[DataDeps.jl](https://github.com/oxinabox/DataDeps.jl):

```julia
using GREBClimate
dir    = greb_data_dir()        # prompts, downloads (~353 MB) and caches on first use
fields = load_greb_jld2!(dir; dataset = :ncep)   # or :era
```

It resolves in this order; only the last step touches the network:

| | Source |
|--:|:-------|
| 1 | an explicit path, `greb_data_dir("/path/to/greb_input_data")` |
| 2 | `ENV["GREB_DATA"]` |
| 3 | `greb_input_data/` beside the package |
| 4 | the `GREB-input-data` DataDep, downloaded and SHA256-verified |

The download is cached under `~/.julia/scratchspaces/.../datadeps/GREB-input-data`,
so it happens once per machine. Two environment variables matter:

| Variable | Effect |
|:---------|:-------|
| `DATADEPS_ALWAYS_ACCEPT=true` | Skip the download prompt. **Required** in CI or any non-interactive session, which otherwise blocks waiting on stdin. |
| `DATADEPS_DISABLE_DOWNLOAD=true` | Make a would-be download throw instead. |

The tests and benchmarks never download: they resolve with
`allow_download=false` and skip data-dependent work when no local copy exists.

## Loading

[`load_greb_jld2!`](@ref) *returns* a [`ClimateFields`](@ref) holding the
climatology, grid geometry, flux corrections and solar table. It does not set
global state: pass the value to [`greb_model!`](@ref) with `fields = ...`. You
can hold several independent instances in one session, e.g. for parameter
sweeps.

A bare `ClimateFields()` is all zeros. Stepping the model on it runs to
completion but produces a meaningless world pinned at the 40 K stability floor
(−233 °C), so `greb_model!` refuses it. Runs that legitimately need no data
opt in explicitly; see [Data-free runs](@ref).

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
│   ├── flux_corrections.jld2       # Tsurf/vapour/Tocean corrections
│   │   # CMIP5 RCP8.5 anomalies - climate-change experiments only
│   ├── cmip5.{tsurf,zonal.wind,meridional.wind,windspeed,omega}.rcp85.ensmean.forcing.jld2
│   │   # ENSO anomalies - :elnino / :lanina only (10 files)
│   └── erainterim.{tsurf,zonal.wind,meridional.wind,windspeed,omega}.{elnino,lanina}.forcing.jld2
├── solar/
│   └── solar_radiation.clim.jld2   # 2D (48×730)
├── solar_scenarios/                # optional
│   ├── solar_paleo.jld2
│   ├── solar_eccentricity.jld2     # (ecc_index, lat, time)
│   └── solar_obliquity.jld2        # (obl_index, lat, time)
└── scenario/                       # optional
    ├── ipcc_scenarios.jld2         # Dict{String,Dict{Int,Float64}}, keyed "rcp85"/"ssp585"/"hist"/...
    └── historical_emissions_population.jld2   # not read by the model (11 KB)
```

## Regenerating the data (maintainers)

Not needed to run the model. The raw GREB `.bin` inputs are collated from
several upstream sources and are not redistributed;
[`DATA_README.md`](https://github.com/EnvDroneSense/GREBClimate.jl/blob/main/DATA_README.md)
documents them.

```bash
julia --project=. tools/convert_greb_to_jld2.jl <input_dir> [output_dir]   # default: greb_input_data/
julia --project=. tools/package_dataset.jl greb_input_data greb_input_data-v1.tar.gz
```

The second command validates the tree against the converter's allowlist,
builds a reproducible archive and prints the SHA256 to paste into
`DATA_SHA256` in `src/data.jl` before attaching the archive to the `data-v1`
release.
