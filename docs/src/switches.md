# Physics Switches

A switch is a [`PhysicsConfig`](@ref) field that turns a process on or off or
selects a parameterization. Build a config with
[`create_experiment_config`](@ref) and override fields afterwards, e.g.
`cfg.log_rain = 1`.

Two suffixes recur:

| Suffix | Meaning | Measures |
|:-------|:--------|:---------|
| `_dmc` | **Deconstruct mean climate**: `false` removes the process from the model | How much the process shapes the mean climate |
| `_drsp` | **Deconstruct response**: `false` keeps the control climate but removes the process's part of the CO₂ response | The process's role in climate sensitivity |

## Overview

| Switch | Type | Default | Family | Effect |
|:-------|:-----|:--------|:-------|:-------|
| `log_clouds_dmc` | Bool | `true` | Mean climate | Zeroes cloud climatology at init when `false` |
| `log_vapor_dmc` | Bool | `true` | Mean climate | No effect - see [Known Limitations](#known-limitations) |
| `log_crcl_dmc` | Bool | `true` | Mean climate | Disables circulation (advection+diffusion+convergence) when `false` |
| `log_hydro_dmc` | Bool | `true` | Mean climate | Disables hydrology (evap/rain/latent heat) when `false` |
| `log_atmos_dmc` | Bool | `true` | Mean climate | Master atmosphere switch - decouples atmosphere from surface when `false` |
| `log_co2_dmc` | Bool | `true` | Mean climate | Forces control-run CO₂ to 0 ppm when `false` |
| `log_ocean_dmc` | Bool | `true` | Mean climate | Master ocean-coupling switch - no mixed-layer/deep-ocean exchange when `false` |
| `log_qflux_dmc` | Bool | `true` | Mean climate | Governs whether flux corrections are computed, loaded, or zeroed |
| `log_clouds_drsp` | Bool | `true` | CO₂ response | Freezes cloud cover at a constant 0.7 when `false` |
| `log_crcl_drsp` | Bool | `true` | CO₂ response | Disables circulation's contribution to the CO₂ response when `false` |
| `log_hydro_drsp` | Bool | `true` | CO₂ response | Disables hydrology's contribution to the CO₂ response when `false` |
| `log_topo_drsp` | Bool | `true` | CO₂ response | Caps topography at 1.0 m ("constant topography") when `false` |
| `log_humid_drsp` | Bool | `true` | CO₂ response | Freezes humidity at a constant 0.0052 kg/kg when `false` |
| `log_ocean_drsp` | Bool | `true` | CO₂ response | Constant mixed-layer depth, no deep-ocean response, when `false` |
| `log_ice` | Bool | `true` | Circulation | Ice-albedo feedback - constant albedo/heat capacity when `false` |
| `log_hdif` | Bool | `true` | Circulation | Horizontal (heat) diffusion |
| `log_hadv` | Bool | `true` | Circulation | Horizontal (heat) advection |
| `log_vdif` | Bool | `true` | Circulation | Vertical/meridional diffusion (water vapor) |
| `log_vadv` | Bool | `true` | Circulation | Advection (water vapor) |
| `log_conv` | Bool | `true` | Circulation | Moisture convergence sub-step (not a `create_experiment_config` keyword) |
| `log_rain` | Int | `0` | Hydrology | Selects the precipitation parameterization |
| `log_eva` | Int | `-1` | Hydrology | Selects the evaporation/exchange-coefficient parameterization |
| `log_clim` | Int | `0` | Hydrology | NCEP-tuned hydrology coefficient override |
| `log_tsurf_ext` | Bool | `false` | External forcing | Loads surface-temperature anomaly forcing |
| `log_hwind_ext` | Bool | `false` | External forcing | Loads horizontal-wind anomaly forcing |
| `log_omega_ext` | Bool | `false` | External forcing | Loads vertical-velocity anomaly forcing |

## Mean-Climate Deconstruction Switches (`_dmc`)

Passed as keywords to `create_experiment_config(:decon_mean_climate; ...)`.

| Switch | Default | Effect when `false` |
|:-------|:--------|:---------------------|
| `log_clouds_dmc` | `true` | Cloud climatology zeroed: no cloud reflection or cloud emissivity. Clouds are prescribed, so this removes them outright. |
| `log_vapor_dmc` | `true` | No effect (dead switch - see [Known Limitations](#known-limitations)). |
| `log_crcl_dmc` | `true` | No transport: no advection, diffusion or convergence. |
| `log_hydro_dmc` | `true` | Humidity climatology zeroed and `hydro!` (evaporation, rain, latent heat) disabled. |
| `log_atmos_dmc` | `true` | Decouples the atmosphere: no sensible heat flux, fixed downward longwave, no transport or hydrology. |
| `log_co2_dmc` | `true` | Control-run CO₂ set to 0 ppm. |
| `log_ocean_dmc` | `true` | No ocean heat storage: land heat capacity everywhere, no sea-ice blending, no deep-ocean exchange. |
| `log_qflux_dmc` | `true` | Together with `log_topo_drsp`, decides whether flux corrections are computed in the spin-up, loaded from file, or zeroed. |

## CO₂-Response Deconstruction Switches (`_drsp`)

Passed as keywords to `create_experiment_config(:decon_2xco2; ...)`, which
doubles CO₂ in the scenario (680 against 340 ppm). There is no `log_co2_drsp`:
the CO₂ doubling is the forcing whose response is being measured.

| Switch | Default | Effect when `false` |
|:-------|:--------|:---------------------|
| `log_clouds_drsp` | `true` | Cloud cover fixed at 0.7 instead of the seasonal climatology. |
| `log_crcl_drsp` | `true` | No transport contribution to the response. |
| `log_hydro_drsp` | `true` | No hydrology contribution to the response. |
| `log_topo_drsp` | `true` | Topography capped at 1 m ("constant topography"). Also affects the flux-correction choice (see `log_qflux_dmc`). |
| `log_humid_drsp` | `true` | Humidity fixed at 0.0052 kg/kg: no water-vapour feedback. |
| `log_ocean_drsp` | `true` | Constant mixed-layer depth and no deep-ocean response: no ocean heat uptake. |

## Circulation Components

Transport uses fixed climatological winds and does not respond to forcing.

| Switch | Default | Effect when `false` |
|:-------|:--------|:---------------------|
| `log_ice` | `true` | No ice-albedo feedback: constant albedo, and heat capacity ignores ice. |
| `log_hdif` | `true` | Disables horizontal (heat) diffusion. |
| `log_hadv` | `true` | Disables horizontal (heat) advection. |
| `log_vdif` | `true` | Disables vertical/meridional diffusion of water vapor. |
| `log_vadv` | `true` | Disables advection of water vapor. |
| `log_conv` | `true` | Disables the moisture-convergence sub-step (water vapor only). Not exposed via `create_experiment_config` - set `cfg.log_conv` directly. |

## Hydrology Parameterization

Evaporation is a bulk formula on the saturation deficit, wind speed and a
prescribed soil wetness. Precipitation scales with humidity, relative humidity
and vertical velocity.

**`log_rain`** - selects the precipitation regression coefficients:

| Value | Meaning |
|:------|:--------|
| `-1` | Original GREB |
| `0` (default) | Best GREB fit (ERA-Interim) |
| `1` | + relative humidity |
| `2` | + omega (vertical velocity) convergence |
| `3` | + relative humidity & omega |

When `log_rain == 1`, a rain-rate limiter is additionally applied.

**`log_eva`** - selects the evaporation/exchange-coefficient parameterization:

| Value | Meaning |
|:------|:--------|
| `-1` (default) | Original GREB (climatological wind + fixed gust term) |
| `0` | `Ts`-derived skin temperature + climatological wind speed + land/ocean exchange coefficients |
| `1` | Modified gust/coefficient variant of `-1` |
| `2` | Modified variant of `0` with different gust/coefficients |

**`log_clim`** - NCEP-vs-ERA coefficient override:

| Value | Meaning |
|:------|:--------|
| `0` (default) | No override |
| `1` | If also `log_rain == 0`, overrides the precipitation coefficients to a hardcoded NCEP-tuned set - independent of which climatology *files* were loaded via `load_greb_jld2!` |

`log_clim = 1` only has an effect together with `log_rain = 0`.

## External-Forcing Gates

Used only by `:elnino`, `:lanina` and `:rcp85`, whose presets already set
all three to `true`. They decide which anomaly fields are loaded from disk.

| Switch | Default | Effect when `true` |
|:-------|:--------|:---------------------|
| `log_tsurf_ext` | `false` | Loads the surface-temperature anomaly forcing field (RCP8.5/ENSO). |
| `log_hwind_ext` | `false` | Loads the horizontal wind (u/v) anomaly field. |
| `log_omega_ext` | `false` | Loads the vertical-velocity (omega) anomaly field. |

## Experiment-Level Fields

Scenario parameters rather than physics switches - set via
`create_experiment_config` keywords or directly on the returned `cfg`:

| Field | Type | Default | Purpose |
|:------|:-----|:--------|:--------|
| `experiment` | `Symbol` | `:full_model` | Selects the experiment branch - see [Experiment Presets](#experiment-presets) below |
| `co2_concentration` | `Float32` | `340.0` | Control-run CO₂ (ppm), and the scenario CO₂ of `:full_model` and `:decon_mean_climate`; `forcing` sets every other experiment's scenario CO₂ |
| `orbital_index` | `Int` | `0` | Row index into the `solar_scenarios` table for `:obliquity`/`:eccentricity` |
| `earth_sun_distance_pct` | `Float64` | `0.0` | Percent change in orbital radius for `:earth_sun_distance` |
| `co2_scenario` | `Dict{Int,Float64}` | `Dict()` | Year→ppm lookup, auto-populated for IPCC RCP/SSP/historical/custom-CO2 experiments |
| `custom_co2_path` | `String` | `""` | User-supplied "year CO2" text file path for `:custom_co2` |

`c_q`, `c_rq`, `c_omega`, and `c_omegastd` also live on `PhysicsConfig` but
are not meant to be set manually - `set_hydrology_parameters!` derives them
from `log_rain`/`log_clim` at the start of every run.

## Experiment Presets

What each `create_experiment_config` preset changes relative to `:full_model`:

| Preset | vs. `:full_model` |
|:-------|:-------------------|
| `:full_model` | Nothing - the baseline |
| `:constant_topo` | `log_topo_drsp = false` |
| `:co2_double` | Scenario at 680 ppm (control stays at 340) |
| `:co2_quadruple` | Scenario at 1360 ppm (control stays at 340) |
| `:solar_plus27` | Solar constant +27 W/m² |
| `:elnino` / `:lanina` | The three `log_*_ext` switches set `true`; adds/subtracts ERA-Interim ENSO anomalies |
| `:paleo_231kyr` | Scenario at 200 ppm with the paleo solar-forcing table (control stays at 340) |
| `:rcp85` | `log_*_ext` switches `true`; loads CMIP5 RCP8.5 anomaly fields |
| `:rcp26` / `:rcp45` / `:rcp60` / `:ssp119` / `:ssp126` / `:ssp245` / `:ssp460` / `:ssp585` | No switches change; loads a year→CO2 lookup table |
| `:historical_co2` | Year counter starts at 1850; loads the observed CO₂ record |
| `:custom_co2` | `custom_co2_path` set from the `co2_path` keyword |
| `:decon_mean_climate` | Exposes `log_clouds_dmc`, `log_ocean_dmc`, `log_atmos_dmc`, `log_co2_dmc`, `log_hydro_dmc`, `log_qflux_dmc`, `log_ice`, `log_hdif`, `log_hadv`, `log_vdif`, `log_vadv` as keywords (all default `true`, i.e. behaves like `:full_model` unless overridden) |
| `:decon_2xco2` | Scenario at 680 ppm (control stays at 340); exposes `log_topo_drsp`, `log_clouds_drsp`, `log_humid_drsp`, `log_ocean_drsp`, `log_hydro_drsp`, `log_ice`, `log_hdif`, `log_hadv`, `log_vdif`, `log_vadv` as keywords |

### Further experiments

`create_experiment_config` covers all of these too. They change no switches
relative to `:full_model` - the experiment symbol alone selects the branch in
`src/tendencies.jl`'s `forcing`, which sets the CO₂ or solar forcing per
timestep. Read that branch for exact behavior.

| Category | Symbols |
|:---------|:--------|
| CO₂ variants | `:a1b_scenario`, `:co2_10x`, `:co2_half`, `:co2_zero`, `:co2_sine_wave`, `:co2_step` |
| Orbital/paleo | `:solar_cycle_11yr`, `:paleo_solar_modern_co2`, `:modern_solar_paleo_co2`, `:obliquity`, `:eccentricity`, `:earth_sun_distance` |
| Other | `:sst_plus1` |
| Regional CO₂ (static mask) | `:regional_co2_nh`, `:regional_co2_sh`, `:regional_co2_tropics`, `:regional_co2_extratropics` |
| Regional CO₂ (dynamic/seasonal mask) | `:regional_co2_ocean`, `:regional_co2_land_ice`, `:regional_co2_winter`, `:regional_co2_summer` |

Two of them take a parameter, passed as a `create_experiment_config` keyword:

| Keyword | Applies to | Meaning |
|:--------|:-----------|:--------|
| `orbital_index` | `:obliquity`, `:eccentricity`, `:paleo_231kyr`, `:paleo_solar_modern_co2` | Which row of the `solar_scenarios/*.jld2` table to load |
| `earth_sun_distance_pct` | `:earth_sun_distance` | Percent change in orbital radius; solar forcing scales as `(1 + 0.01·pct)^-2` |

The 16 `log_*` keywords apply **only** to `:decon_mean_climate` and
`:decon_2xco2`. Passing one for any other experiment emits a warning and is
ignored - set the field on the returned `PhysicsConfig` instead.

`:regional_co2_ocean` and `:regional_co2_land_ice` derive their mask from the
control run's annual-mean ice cover, so it is built once per run by
`apply_dynamic_co2_mask!` between the control and scenario phases rather than
by `forcing`. The control run itself always sees an all-ones mask.

## Known Limitations

- **`log_vapor_dmc`** is declared on `PhysicsConfig` but is never read
  anywhere in `src/` - setting it currently has no effect on the model.
- **`log_conv`** works like the other circulation-component switches
  (gates the moisture-convergence sub-step) but, unlike them, is not
  exposed as a keyword by `create_experiment_config` - it must be set with
  `cfg.log_conv = ...` after construction.

## References

1. Dommenget, D., and Flöter, J. (2011). Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model. *Climate Dynamics*, 37: 2143. [doi:10.1007/s00382-011-1026-0](https://doi.org/10.1007/s00382-011-1026-0)
2. Stassen, C., Dommenget, D., and Loveday, N. (2019). A hydrological cycle model for the Globally Resolved Energy Balance (GREB) model v1.0. *Geoscientific Model Development*, 12, 425-440. [doi:10.5194/gmd-12-425-2019](https://doi.org/10.5194/gmd-12-425-2019)

See the repository [README](https://github.com/EnvDroneSense/GREBClimate.jl#references)'s References section for the original GREB model homepage link.
