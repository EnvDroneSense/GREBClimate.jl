# Configuration

A run is described by a [`Config`](@ref) and its length by a
[`RunSpec`](@ref). Build a `Config` from a named experiment with
[`preset`](@ref), or from its parts:

| Part | Type | Decides |
|:-----|:-----|:--------|
| `scenario` | [`Scenario`](@ref) | What the experiment imposes: CO₂ path and where it applies, sunlight, surface forcing, control CO₂, start year, output form |
| `processes` | [`Processes`](@ref) | Which physical processes run, each with a stated meaning when off |
| `hydrology` | [`Hydrology`](@ref) | The evaporation and rain scheme |
| `corrections` | [`Corrections`](@ref) | How the flux corrections are obtained: [`SpinUp`](@ref)`(years)`, [`Stored`](@ref)`()`, [`NoCorrections`](@ref)`()` |

```julia
config = preset(:co2_double)                                   # a named experiment
config = preset(:co2_double; processes = (clouds = :uniform,))  # with one process changed
config = Config(scenario = Scenario(co2 = ConstantCO2(500), solar = SolarConstant(10)))
```

Every part validates its options when it is built, so a misspelt option fails
at once instead of being ignored.

## Processes

| Option | Values (default first) | Meaning of the other values |
|:-------|:-----------------------|:----------------------------|
| `atmosphere` | `true`, `false` | No longwave back-radiation, sensible or latent heat, or transport |
| `clouds` | `:observed`, `:none`, `:uniform` | Cloud cover 0, or 0.7 everywhere |
| `humidity` | `:observed`, `:uniform` | Humidity climatology 0.0052 kg/kg everywhere |
| `hydrology` | `:full`, `:no_evap_rain`, `:none` | No evaporation and rain; or no water cycle (humidity 0, not updated) |
| `ocean` | `:full`, `:mixed_layer`, `:none` | A 50 m mixed layer without deep ocean; or land heat capacity everywhere |
| `topography` | `:observed`, `:flat` | Topography capped at 1 m |
| `co2` | `true`, `false` | CO₂ 0 ppm in the control and the scenario |
| `ice_albedo` | `true`, `false` | No ice-albedo feedback |
| `transport` | `true`, `false` | No atmospheric heat and moisture transport |
| `heat_diffusion`, `heat_advection`, `vapour_diffusion`, `vapour_advection`, `moisture_convergence` | `true`, `false` | That transport term off |

## Hydrology

| Option | Values (default first) |
|:-------|:-----------------------|
| `rain` | `:fitted` (Stassen et al. 2019), `:original` (proportional to humidity), `:rh` (+ relative humidity, with a rain-rate limit), `:omega` (+ vertical velocity), `:rh_omega` |
| `rain_fit` | `:era`, `:ncep`: the reanalysis the `:fitted` coefficients were fitted to; independent of the climatology files [`load_greb_jld2!`](@ref) reads |
| `evaporation` | `:original` (climatological wind plus a fixed gust term), `:skin` (skin temperature, land/ocean exchange coefficients), `:original_gust`, `:skin_gust` (the same with modified gust terms) |

[`mscm_hydrology`](@ref)`()` is the original GREB scheme; with
`Processes(moisture_convergence = false)` it reproduces the MSCM model.

## Presets

[`preset_names`](@ref)`()` lists them. Every preset runs the full model with a
3-year spin-up unless stated.

| Preset | Scenario |
|:-------|:---------|
| `:full_model` | 340 ppm throughout |
| `:co2_double`, `:co2_quadruple`, `:co2_10x`, `:co2_half`, `:co2_zero` | 680, 1360, 3400, 170, 0 ppm |
| `:co2_sine_wave`, `:co2_abrupt_reverse`, `:a1b` | [`CO2SineWave`](@ref), [`CO2Step`](@ref)`(680, 340, 1980)`, [`A1BRamp`](@ref) (control at 280 ppm) |
| `:rcp26`, `:rcp45`, `:rcp60`, `:rcp85`, `:ssp119`, `:ssp126`, `:ssp245`, `:ssp460`, `:ssp585`, `:historical_co2` | CO₂ from the dataset's table (control at 280 ppm; `:historical_co2` starts in 1850) |
| `:custom_co2` | CO₂ from a `year ppm` file: `preset(:custom_co2; path = "co2.txt")` |
| `:solar_plus27`, `:solar_cycle_11yr` | Solar constant +27 W/m², or an 11-year cycle of ±1 W/m² |
| `:paleo_231kyr`, `:paleo_solar_modern_co2`, `:modern_solar_paleo_co2` | Insolation 231 kyr ago and/or 200 ppm |
| `:obliquity`, `:eccentricity`, `:earth_sun_distance` | Orbital insolation tables (`index = ...`, default the row nearest today) or a distance change (`pct = ...`); output absolute, start year 1 |
| `:elnino`, `:lanina`, `:rcp85_boundary` | Boundary anomalies of surface temperature, winds and vertical velocity at 340 ppm (`:rcp85_boundary`: the CMIP5 RCP8.5 change; `:rcp85` is the CO₂ path) |
| `:sst_plus1` | Ocean surface held at climatology + 1 K |
| `:regional_co2_nh`, `_sh`, `_tropics`, `_extratropics`, `_ocean`, `_land_ice`, `_winter`, `_summer` | 680 ppm in one region or season, 340 ppm elsewhere |
| `:decon_mean_climate` | Mean-climate deconstruction: 340 ppm, MSCM physics, stored corrections. Switch processes off with `processes = (...)`; run with `RunSpec(scnr = 0)` |
| `:decon_2xco2` | 2×CO₂-response deconstruction: 680 ppm, default physics |

On the stored corrections a switched-off process changes the climate; a
spin-up would recompute the corrections for each configuration and pull every
one back to the observed climate.

## References

1. Dommenget, D., and Flöter, J. (2011). Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model. *Climate Dynamics*, 37: 2143. [doi:10.1007/s00382-011-1026-0](https://doi.org/10.1007/s00382-011-1026-0)
2. Stassen, C., Dommenget, D., and Loveday, N. (2019). A hydrological cycle model for the Globally Resolved Energy Balance (GREB) model v1.0. *Geoscientific Model Development*, 12, 425-440. [doi:10.5194/gmd-12-425-2019](https://doi.org/10.5194/gmd-12-425-2019)

See the repository [README](https://github.com/EnvDroneSense/GREBClimate.jl#references)'s References section for the original GREB model homepage link.
