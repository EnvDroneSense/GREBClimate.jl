# Model overview

GREB is a conceptual climate model: it solves energy and moisture budgets on a
fixed 3.75° × 3.75° grid (96 × 48), but takes the atmospheric circulation,
clouds, soil moisture and ocean mixed-layer depth from observed seasonal
climatologies instead of simulating them; the circulation does not respond to
a forcing. That keeps it fast enough for many simulated years while still
resolving where a forcing acts (Dommenget & Flöter, 2011).

## Prognostic variables

| Variable | Meaning | Unit |
|:---------|:--------|:-----|
| `Ts` | Surface temperature (land, or ocean mixed layer) | K |
| `Ta` | Atmospheric temperature | K |
| `To` | Deep-ocean temperature, the slow reservoir | K |
| `q` | Specific humidity of the atmosphere | kg/kg |

## Energy balance

The surface temperature changes with the sum of the heat fluxes into it:

```math
c_{surf} \frac{dT_{surf}}{dt} = F_{solar} + F_{thermal} + F_{latent} + F_{sense} + F_{ocean} + F_{correct}
```

- `c_surf`: surface heat capacity - the ocean mixed layer over ice-free ocean,
  2 m of soil over land, falling to 2 m of water where sea ice forms
- `F_solar`, `F_thermal`: absorbed sunlight, and longwave emission minus back-radiation
- `F_latent`, `F_sense`: evaporative cooling, and turbulent heat exchange with the atmosphere
- `F_ocean`: heat exchange with the deep ocean
- `F_correct`: flux correction (see [A run](@ref))

The atmospheric temperature and humidity have similar budgets, plus transport
by the circulation.

## Components

Each timestep (12 hours, 730 per year) [`time_loop!`](@ref) calls
[`tendencies!`](@ref) for the radiation, hydrology, circulation and deep-ocean
terms, integrates the four prognostic fields, then updates sea ice.

| Component | Function | What it computes |
|:----------|:---------|:-----------------|
| Shortwave radiation | [`SWradiation!`](@ref) | Ice cover and surface albedo from `Ts` (albedo rises linearly as the surface cools through a band just below freezing), then the absorbed solar flux. Cloud albedo scales with the ISCCP cloud-cover climatology. |
| Longwave radiation | [`LWradiation!`](@ref) | Atmospheric emissivity from CO₂, water vapour and clouds (a 10-parameter log fit), then the up- and downward longwave fluxes. This is where the greenhouse effect lives. |
| Hydrology | [`hydro!`](@ref) | Evaporation by a bulk formula (four variants, `log_eva`), precipitation from humidity, relative humidity and vertical velocity (`log_rain`), and the latent heat flux (Stassen et al., 2019). |
| Circulation | [`circulation!`](@ref) | Transport of `Ta` and `q` by advection with the climatological 850 hPa winds and isotropic diffusion, plus moisture convergence for `q` from the vertical-velocity climatology. Runs 24 half-hour sub-steps per timestep, more near the poles; about 93% of the run time. |
| Deep ocean | [`deep_ocean!`](@ref) | Heat exchange between the mixed layer (`Ts`) and the deep ocean (`To`) by entrainment, detrainment and turbulent mixing. |
| Sea ice | [`seaice!`](@ref) | The surface heat capacity where sea ice forms; latent heat of freezing is neglected. |
| Forcing | [`forcing`](@ref) | The scenario's CO₂ and solar multiplier for the current timestep, per experiment. |

## Shortwave add-ons

The solar multiplier is one value per latitude row. After [`forcing`](@ref)
sets it, two optional add-ons on the config multiply into it, in the scenario
run only: a solar record (`cfg.solar_scenario`, from
[`load_solar_series`](@ref)) and a stratospheric aerosol
(`cfg.aerosol`, an [`AerosolScenario`](@ref)). Anomalies are therefore relative
to an undimmed control. Both are off by default and leave default runs
bit-identical.

The aerosol path runs eruptions, sustained injections or a published record
through three steps:

| Step | Method | Main assumptions |
|:-----|:-------|:-----------------|
| Optical depth in time | Two linear reservoirs: peak 5 months after injection, e-folding decay of 1 year | Tropical Pinatubo values, used for every injection class |
| Optical depth in latitude | Fixed profile per class: tropical, northern or southern extratropical | Tropical spreads into both hemispheres; high latitudes keep 80% of the tropical value |
| Shortwave transmission | Delta-Eddington layer (Joseph et al., 1976) at a fixed zenith angle, on a calibrated optical depth | Calibrated to −30 W/m² per unit optical depth (Sato et al., 1993); uncalibrated, the layer gives −82 W/m² |

The layer is a multiplier on incoming sunlight: it does not reflect light
coming back up from clouds and surface, which is the main reason for the
calibration. Not represented: the aerosol's longwave effect, stratospheric
heating, ozone chemistry, particle growth over an eruption, and poleward
transport, so mid-latitudes peak at the same time as the tropics instead of
3 to 8 months later. The parameters are provisional.

Compared with observed temperatures (GISTEMP, with trend, ENSO and the annual
cycle removed), the model cools about twice as strongly as observed after
Agung (1963), El Chichón (1982) and Pinatubo (1991), and recovers more slowly.
The ratio is similar for all three, so comparisons between aerosol scenarios
are meaningful, but the absolute cooling is too large. Whether the excess lies
in the forcing or in the model's response to it is not yet known. The
comparison script is `tools/validate_volcanic.jl`.

## A run

[`greb_model!`](@ref) runs up to three phases, set by [`RunSpec`](@ref):

| Phase | CO₂ | Purpose |
|:------|:----|:--------|
| Flux-correction spin-up (`flux` years) | control | [`qflux_correction!`](@ref) derives the corrections for `Ts`, `To` and `q` that hold the control at the observed climatology; without them the model drifts by several K |
| Control (`ctrl` years) | control: 340 ppm (280 for the IPCC scenarios) | The reference climate |
| Scenario (`scnr` years) | set by [`forcing`](@ref) per experiment | The experiment |

The result holds monthly means (`MonthlyRecord`s) of 13 fields for the control
and the scenario. The scenario is returned as an anomaly: each month minus the
same calendar month of the control's final year, except for the orbital
experiments and runs without a control. See the [Tutorial](@ref) for a run
and the [Physics Switches](@ref) for the switches each component reads.

## References

- Dommenget & Flöter (2011), *Climate Dynamics* 37: 2143 - the model and its
  energy balance.
- Stassen, Dommenget & Loveday (2019), *Geosci. Model Dev.* 12: 425 - the
  hydrological cycle.
- Joseph, Wiscombe & Weinman (1976), *J. Atmos. Sci.* 33: 2452 - the
  delta-Eddington approximation.
- Sato, Hansen, McCormick & Pollack (1993), *J. Geophys. Res.* 98: 22987 -
  stratospheric aerosol optical depths and their forcing.
