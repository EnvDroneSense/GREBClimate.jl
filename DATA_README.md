# GREB raw input data (maintainers)

> **You do not need this to run the model.** GREBClimate reads a prepared
> `.jld2` dataset that it downloads on first use; see
> [Input data](https://EnvDroneSense.github.io/GREBClimate.jl/dev/data/).
> This page lists the **raw** GREB `.bin` files that dataset is generated from.
> They are collated from several upstream sources (NCEP, ERA-Interim, ISCCP,
> WOCE, CMIP5, IPCC scenario tables) and are not redistributed here.

## Converting

```bash
julia --project=. tools/convert_greb_to_jld2.jl <input_dir> [output_dir]   # defaults: Data/, greb_input_data/
julia --project=. tools/package_dataset.jl greb_input_data greb_input_data-v1.tar.gz
```

## Generating new input data (ERA5)

`tools/fetch_era5_data.py` builds an ERA5 climatology from the Copernicus
Climate Data Store and writes it in the `.bin` format the converter reads.
It requires your own [CDS API access](https://cds.climate.copernicus.eu/how-to-api).

```bash
python -m pip install cdsapi netCDF4 numpy
python tools/fetch_era5_data.py --start-year 1991 --end-year 2020 --out-dir era5_raw
julia --project=. tools/convert_greb_to_jld2.jl era5_raw greb_input_data_era5

The converter reads the `.bin` files (and matching `.ctl` files, if present)
**flat** from the input directory:

```
Data/                                  # <input_dir>
├── *.bin, *.ctl                       # the fields below
├── ipcc.scenario.*.forcing*.txt       # CO2 scenario tables (optional)
└── solar_forcing_scenarios/           # orbital/paleo solar tables (optional)
    └── greb.solar.*.bin
```

It converts only the fields the model reads - `MODEL_FIELD_NAMES` in
`tools/convert_greb_to_jld2.jl` is the authoritative list - and warns about any
that are missing. `--all` converts every `.bin` present. `package_dataset.jl`
checks the result against the same list before building the release archive.

## Required fields

All 3D fields are 96×48×730 (one year at 12-hour steps).

**Static (96×48)**

| File | Content |
|------|---------|
| `global.topography.bin` | Topography (m); ocean points < 0 |
| `greb.glaciers.bin` | Glacier mask (0 or 1) |

**Solar (48×730)**

| File | Content |
|------|---------|
| `solar_radiation.clim.bin` | 24-hour mean top-of-atmosphere solar radiation (W/m²) |

**Climatology: one of the two datasets** (`dataset=:ncep` or `:era`)

| NCEP | ERA-Interim | Content |
|------|-------------|---------|
| `ncep.tsurf.1948-2007.clim.bin` | `erainterim.tsurf.1979-2015.clim.bin` | Surface temperature (K) |
| `ncep.zonal_wind.850hpa.clim.bin` | `erainterim.zonal_wind.850hpa.clim.bin` | Zonal wind at 850 hPa (m/s) |
| `ncep.meridional_wind.850hpa.clim.bin` | `erainterim.meridional_wind.850hpa.clim.bin` | Meridional wind at 850 hPa (m/s) |
| `ncep.atmospheric_humidity.clim.bin` | `erainterim.atmospheric_humidity.clim.bin` | Specific humidity (kg/kg) |
| `ncep.soil_moisture.clim.bin` | (uses the NCEP file) | Soil moisture fraction (0-1) |

**Climatology: used with both datasets**

| File | Content |
|------|---------|
| `isccp.cloud_cover.clim.bin` | Cloud cover fraction (0-1) |
| `woce.ocean_mixed_layer_depth.clim.bin` | Ocean mixed-layer depth (m) |
| `Tocean.clim.bin` | Deep-ocean temperature (K) |
| `erainterim.omega.vertmean.clim.bin` | Vertical velocity (Pa/s) |
| `erainterim.omega_std.vertmean.clim.bin` | Standard deviation of vertical velocity |
| `erainterim.windspeed.850hpa.clim.bin` | Wind speed at 850 hPa (m/s) |

**Flux corrections** (combined into one `flux_corrections.jld2`)

| File | Content |
|------|---------|
| `Tsurf_flux_correction.bin` | Surface temperature correction (W/m²) |
| `vapour_flux_correction.bin` | Water vapour correction (kg/m²/s) |
| `Tocean_flux_correction.bin` | Deep-ocean correction (W/m²) |

**Anomaly forcing** (only for the experiments named)

| Files | Used by |
|-------|---------|
| `cmip5.{tsurf,zonal.wind,meridional.wind,windspeed,omega}.rcp85.ensmean.forcing.bin` | `:rcp85` (CMIP5 ensemble-mean anomalies) |
| `erainterim.{tsurf,zonal.wind,meridional.wind,windspeed,omega}.{elnino,lanina}.forcing.bin` | `:elnino`, `:lanina` |

## Optional inputs

**CO2 scenario tables.** Whitespace-separated text, one row per year, no
header: `year CO2`. Extra columns are ignored. They are combined into
`scenario/ipcc_scenarios.jld2`, keyed by the name between `ipcc.scenario.`
and `.forcing` (e.g. `"rcp85"`, `"hist"`).

| File | Scenario |
|------|----------|
| `ipcc.scenario.rcp{26,45,6,85}.forcing.txt` | RCP 2.6, 4.5, 6.0, 8.5 |
| `ipcc.scenario.ssp{119,126,245,460,585}.forcing.txt` | SSP1-1.9, 1-2.6, 2-4.5, 4-6.0, 5-8.5 |
| `ipcc.scenario.hist.forcing.CO2.emission.pop.txt` | Historical CO2, 1850-2017 |

The CO2 values in the RCP files are GREB's simplified forcing index, not
literal atmospheric ppm. `:historical_co2` starts its clock at 1850 (the other
scenarios at 1950), so `RunSpec(scnr=168)` covers the whole record. The
historical file's columns 3-4 (emissions, population) are not read by the
model; the converter keeps them in `scenario/historical_emissions_population.jld2`
for other use.

**Solar forcing scenarios** in `solar_forcing_scenarios/`:

| Files | Experiment |
|-------|------------|
| `greb.solar.eccentricity.{0-60}.bin` | `:eccentricity` |
| `greb.solar.obliquity.{0-230}.bin` (steps of 5) | `:obliquity` |
| `greb.solar.231K_hybers.corrected.bin` | `:paleo_231kyr` |

## File format

GrADS binary: 32-bit little-endian floats, no header, Fortran (column-major)
order longitude × latitude × time. Longitude and latitude are 96 × 48 points at
3.75°. The `.ctl` files describe this layout for GrADS; the converter copies
their text into the dataset but does not need them.

## Troubleshooting

| Problem | Check |
|---------|-------|
| The converter warns about missing fields | File names match exactly (case-sensitive on Linux/macOS) and sit flat in the input directory |
| A field fails to convert | Its size is 96×48, 48×730 or 96×48×730 Float32 values with no header bytes |
| Unrealistic model output | The flux-correction files are present, and topography is negative over ocean |
