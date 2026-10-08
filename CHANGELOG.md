# Changelog

Notable changes to GREBClimate.jl, following
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Changed

- **Breaking:** the fields of `MonthlyAccumulator` take the names of
  `MonthlyRecord` (`Tmm` is now `Ts`, `apmm` is `albedo`, `olrmm` is `olr`, and
  so on), and its keyword constructor is gone; `MonthlyAccumulator()` is
  unchanged. The output fields are declared once, as
  `GREBClimate.OUTPUT_FIELDS`. Results are unchanged.
- **Breaking:** one time value, `ModelTime(phase, start_year, step)`, says
  where a run is; day, month and year are derived from it. Results are
  unchanged. What changes for a caller:
  - `greb_model!` returns `(ctrl, scnr, ctrl_time, scnr_time, scnr_anomaly)`.
    `ctrl_time` and `scnr_time` hold the `(year, month)` of each record;
    `scnr_anomaly` is `true` when `scnr` is the anomaly against the control and
    `false` when it holds absolute values.
  - `forcing(t::ModelTime, r)` replaces `forcing(it, year, r)`.
  - `scenario_anomalies(scnr, scnr_time, ctrl_clim)` takes the record times
    and no longer assumes that the first record is January.
  - The observer's view has `time`, `phase`, `step`, `year` and
    `step_of_year` in place of `phase`, `it`, `year` and `ityr`. `phase` is
    `GREBClimate.control` or `GREBClimate.scenario`, no longer `:ctrl` or
    `:scnr`; `RangeCheck`'s `first` has `step` in place of `it`.
  - `time_loop!`, `output!` and `diagnostics!` take the time value and return
    `nothing`; the month and record counters are gone.
- **Breaking:** `TimeState` is removed. The physics kernels (`SWradiation!`,
  `LWradiation!`, `hydro!`, `seaice!`, `deep_ocean!`, `convergence!`,
  `advection!`, `circulation!`) and `tendencies!` take the index of the
  climatology slice, an `Int`, where they took a `TimeState`; `time_loop!` and
  `qflux_correction!` lose the argument. Results are unchanged.
- **Breaking:** `ClimateFields` no longer carries the ten anomaly arrays
  (`Ts_clim_anom_enso`, `Ts_clim_anom_cc` and the others). The anomalies of a
  `BoundaryAnomaly` scenario are in `fields.boundary_anomaly`, which is
  `nothing` until `load_boundary_anomaly!` loads a source and then holds five
  arrays for one source at a time; `anom_cc_source` and `anom_enso_source`
  become `boundary_anomaly_source`. A `ClimateFields` is 135 MB smaller, also
  per ensemble member. Results are unchanged.
- Internal calendar constants are written out: `ndt_days`, `ndays_yr`,
  `cjday_mon` and `jday_mon_cumsum` are `steps_per_day`, `days_per_year`,
  `days_in_month` and `last_day_of_month`.
- The plotting toolbox (`viz/`) reads the year and month of each record from
  the result's `ctrl_time` and `scnr_time`, and whether the scenario is an
  anomaly from its `scnr_anomaly` (no longer guessed from the temperature):
  animation frames are titled
  `Jul 1950`, and yearly means and the seasonal cycle group by the record
  times. `annual`, `seasonal_cycle` and `map_frames` take the times as an
  optional argument; a bare record vector is counted from January of year 1,
  as before.

### Fixed

- `:a1b` (`A1BRamp`) holds 700 ppm after 2100, as the Fortran does. It fell
  back to 340 ppm in 2101, so a run of more than 151 scenario years lost its
  forcing. Runs that end by 2100 are unchanged.

## [2.1.0] - 2026-10-05

### Added

- Calendar functions `step_of_year`, `day_of_year`, `month_of_step` and
  `decimal_year`. The calendar is computed in one place and the 200-year
  lookup table is gone (2.3 MB less memory); results are unchanged.
- `run_ensemble([reduce,] run, configs; fields, jld2_dir, ntasks, logger)`: runs
  many configurations side by side as tasks, each on its own copy of the
  fields, and returns one entry per member. Every member gives exactly the
  result of running it alone.
- `global_mean(field)`: the area-weighted global mean of a field on the model
  grid. The examples use it; their plain means read 8 to 11 K too cold.
- `GREBClimate.RangeCheck()`: an observer that records the first step at which
  `Ts`, `Ta`, `To` or `q` leaves a physical range, for a run that diverges
  without going non-finite.
- `load_climatology(dir; corrections = false)` loads a dataset without the
  stored flux corrections; the correction arrays stay zero.

### Changed

- **Breaking:** the loaders no longer carry the file format in their names,
  and the two anomaly loaders are one function that takes the source:

  | Before | Now |
  |---|---|
  | `load_greb_jld2!(dir)` | `load_climatology(dir)` |
  | `load_flux_corrections_jld2!(dir, fields)` | `load_flux_corrections!(dir, fields)` |
  | `load_co2_scenario_jld2(dir, key)` | `load_co2_scenario(dir, key)` |
  | `load_custom_co2_scenario(path)` | `load_co2_custom(path)` |
  | `load_solar_forcing_jld2(dir, kind, index)` | `load_solar_forcing(dir, kind, index)` |
  | `load_cc_anomaly_jld2!(dir, fields)` | `load_boundary_anomaly!(dir, fields, :cmip5_rcp85)` |
  | `load_enso_anomaly_jld2!(dir, fields, event)` | `load_boundary_anomaly!(dir, fields, event)` |
  | `read_jld2(path)` | `read_field(path)` |

- **Breaking:** the `Processes` switches `vapour_diffusion` and
  `vapour_advection` are `vapor_diffusion` and `vapor_advection`, the
  spelling of every other name in the package.
- **Breaking:** the fields of `CirculationWorkspace` have no `_buf` or `_out`
  suffix: `ws.Q_lat_buf` is `ws.Q_lat`, `ws.precip_out` is `ws.precip`, and
  so on for all 24. The type itself is `ModelWorkspace`; it was
  `CirculationWorkspace`, although every kernel uses it.
- **Breaking:** the fields of `ClimateFields` separate their words:

  | Before | Now |
  |---|---|
  | `Tclim`, `Toclim`, `qclim` | `Ts_clim`, `To_clim`, `q_clim` |
  | `uclim`, `vclim`, `wsclim` | `u_clim`, `v_clim`, `wind_speed_clim` |
  | `omegaclim`, `omegastdclim` | `omega_clim`, `omega_std_clim` |
  | `mldclim`, `cldclim`, `swetclim` | `mld_clim`, `cloud_clim`, `soil_wetness_clim` |
  | `TF_correct`, `qF_correct`, `ToF_correct` | `Ts_flux_correction`, `q_flux_correction`, `To_flux_correction` |

  The anomaly and split-wind fields follow their stem (`Ts_clim_anom_cc`,
  `u_clim_pos`). `ModelState.Tsmn` is `Ts_annual_mean`, and
  `ResolvedHydrology.c_omegastd` is `c_omega_std`.
- **Breaking:** scenario fields say what they hold: `SSTOffset.offset` (was
  `K`), `SolarConstant.offset` (was `dW`), `EarthSunDistance.percent` (was
  `pct`), `SeasonalCO2.in_season` and `.out_of_season` (were `inside` and
  `outside`). The keyword of the preset follows:
  `preset(:earth_sun_distance; percent = 1.5)`. Positional construction is
  unchanged.
- **Breaking:** `build_monthly_climatology` is `monthly_climatology`,
  `apply_scenario_anomalies` is `scenario_anomalies`,
  `compute_annual_ice_climatology` is `ice_climatology`, and
  `apply_dynamic_co2_mask!` is `apply_surface_mask!`.
- **Breaking:** the kernels and loop functions are no longer exported, only
  what a user calls is. `SWradiation!`, `LWradiation!`, `hydro!`,
  `convergence!`, `seaice!`, `deep_ocean!`, `diffusion!`, `advection!`,
  `circulation!`, `tendencies!`, `time_loop!`, `output!`, `diagnostics!`,
  `qflux_correction!` and `init_model!` are reached as `GREBClimate.name`
  or with `using GREBClimate: name`. Their names and behaviour are unchanged.
- **Breaking:** arguments that were never read are gone. `diagnostics!` is
  `diagnostics!(year, surf, state, timestate)` and `diffusion!` is
  `diffusion!(T1, h_scl, fields, ws)`.
- **Breaking:** `ModelWorkspace` and `MonthlyAccumulator` are immutable. Their
  arrays are written in place as before; a field can no longer be replaced.
- Internal constants that take part in the model's formulas are lower case
  (`is_polar`, `polar_diff_time2`, `ΔT_air_factor`); capitals are kept for
  the package's own settings and tables. Some are named after what they are:
  `pi_f32` (was `const_pi`), `convergence_factor` (`const_factor`),
  `q_to_mm_per_day` (`conv_factor`), `solar_percent` (`S0_var`) and
  `emissivity_fit` (`p_emi`). No exported name is affected.
- What a run and the loaders report goes through the logger (`@info`) and no
  longer through `println`: the lines carry an `[ Info:` prefix, go to
  standard error, and are silenced with
  `with_logger(NullLogger())` instead of `redirect_stdout(devnull)`.
- The precompile workload runs a scenario year as well as a control year, so
  the first scenario run of a session compiles less.
- The yearly line names what it shows:
  `1970: Ts global mean 13.82 °C; 178 E 9 N 26.79; 58 E 51 N 4.85`. A spin-up
  year reads `spin-up year 1: ...`, and a run without a scenario no longer
  reports a scenario of 0 years.
- **Breaking:** the RCP6.0 CO2 table is `CO2Table(:rcp60)`, the name of its
  preset. `load_co2_scenario(dir, :rcp6)`, and so a `CO2Table(:rcp6)` run, raise an
  `ArgumentError`. The `:rcp60` preset is unchanged.
- `load_climatology` and `load_flux_corrections!` raise an `ArgumentError`
  when `climatology/flux_corrections.jld2` or one of its three tables is
  missing. They used to fill the corrections with zeros and warn. For a
  dataset without the file, pass `corrections = false`.
- A `BoundaryAnomaly` scenario (`:rcp85_boundary`, `:elnino`, `:lanina`) reads
  its anomaly files once per `ClimateFields`: a later run on the same fields,
  directory and source reuses the arrays. `ClimateFields` has two new fields
  for this, `anom_cc_source` and `anom_enso_source`. Results are unchanged.

## [2.0.1] - 2026-10-04

### Changed

- A scenario with an `SSTOffset` (`:sst_plus1`) no longer allocates a copy of
  the temperature climatology at every step. Results are unchanged.
- The `observer` keyword of `greb_model!` and `GREBClimate.BudgetCheck` are no
  longer marked experimental. The contents of the observer's `view` are part
  of the interface from here on.

### Changes to model results

- `:regional_co2_nh`, `:regional_co2_sh`, `:regional_co2_tropics` and
  `:regional_co2_extratropics` (any scenario with a `LatitudeMask`): the mask
  was already active in the spin-up and the control, which therefore ran on
  half the control CO2 outside the band, and the scenario then doubled the CO2
  everywhere relative to that control. The mask now applies from the start of
  the scenario only, as in the original GREB code (`subroutine forcing`). **The
  results of these four presets change**: the response outside the band is
  much smaller than before. `SurfaceMask` presets are unaffected.

### Removed

- Seven fields of `CirculationWorkspace` that were written every step and read
  by nothing: `qs`, `rq`, `Tskin`, `ws_base`, `cE_buf`, `temp_buf` and
  `a_atmos_buf`. Results are bit-identical.

### Fixed

- `corrections = Stored()` without `jld2_dir` no longer zeroes the flux
  corrections that `load_greb_jld2!` put in `fields`: it runs on them. Before,
  such a run (for example `greb_model!(run, preset(:decon_mean_climate);
  fields)`) logged one warning and ran without corrections. With `jld2_dir`
  given, a missing corrections file is now an `ArgumentError`.
- `load_greb_jld2!` with an unknown `dataset` (for example `:era5`) raises an
  `ArgumentError` naming the valid ones. Before, it loaded the NCEP files and
  printed the name it was given.
- A scenario with a `SurfaceMask` (`:regional_co2_ocean`,
  `:regional_co2_land_ice`) run with `RunSpec(ctrl = 0)` raises an
  `ArgumentError`. Before, the mask was built from an ice cover of zero.
- Several `greb_model!` runs in parallel tasks no longer fail with a JLD2
  error (`InvalidDataException`, `EOFError`) when they read dataset files at
  the same time: the loaders take a lock around each file open.
- `load_solar_forcing_jld2` with an orbital `index` that is not in the table
  raises an `ArgumentError` listing the available ones, and `load_greb_jld2!`
  raises an error naming the file when the solar table has the wrong shape.
  Both were `@assert`s, which Julia may skip.
- Docstrings corrected: `CO2Table` lists the `:rcp85` key; `Hydrology`
  describes `:skin_gust` as it is computed (no skin temperature);
  `ClimateFields` no longer says one instance per run, and states that
  all-zero fields give NaN output, not a 40 K world; `tendencies!` and
  `mscm_hydrology` had a cut-off sentence and a missing separator.

## [2.0.0] - 2026-10-01

### Added

- A new configuration API: `Config` with a `Scenario` (CO2 path and mask,
  sunlight, surface forcing, control CO2, start year, output form), a
  `Processes` set (one option per process, each with an explicit meaning when
  off, replacing the overlapping mean-climate and 2xCO2 switch sets), a
  `Hydrology` scheme (named options instead of integer codes; `mscm_hydrology()`
  for the original GREB scheme) and `Corrections` (`SpinUp(years)`, `Stored()`,
  `NoCorrections()`). `preset(name)` builds every experiment; run it with
  `greb_model!(run, config)`, where the spin-up length comes from `SpinUp`.
  Results are identical to the former experiments. Any topography works with
  any corrections; every preset spins up for 3 years unless told otherwise,
  except `:decon_mean_climate` (below). `preset(name; processes = (ocean =
  :none,))` changes single options of a preset's physics. Scenario parts
  combine freely (for example `Scenario(co2 = ConstantCO2(500), solar =
  SolarConstant(10))`). `resolve(config)` returns a `ResolvedConfig` with the
  CO2 or solar table and the rain coefficients the configuration refers to.
  Preset names that changed: `:rcp85` is `:rcp85_boundary`, `:co2_step` is
  `:co2_abrupt_reverse`, `:a1b_scenario` is `:a1b`; `:constant_topo` has no
  preset (use `preset(:co2_double; processes = (topography = :flat,),
  corrections = Stored())`).
- `:rcp85`: RCP8.5 CO2 from the dataset's table, like the other IPCC presets
  (control at 280 ppm). The boundary-anomaly run is `:rcp85_boundary`.
- A regression test that the MSCM configuration (`mscm_hydrology()`,
  `moisture_convergence = false`) reproduces the MSCM 2xCO2 response: year 1
  global mean 0.5946 K in both.
- `MonthlyRecord` has two more fields: `olr`, the longwave leaving to space
  (positive upward), and `lwdown`, the longwave the air sends to the surface
  (positive downward), both in W/m2. Code that builds a `MonthlyRecord` by hand
  or relies on its 13 fields has to add them.
- `greb_model!(...; observer = f)` (experimental): `f(point, view)` is called
  before and after every step's state update of the control and scenario runs
  with the model's state and flows, for diagnostics that need per-step values.
  Without it the run is unchanged. `GREBClimate.BudgetCheck()` is an observer
  that checks each store changes by the sum of its flows and counts the cells
  the model's limiters held; `tools/diagnostics/budget.jl` prints its result
  and the global-mean flows of a run.

### Removed

- `PhysicsConfig`, `create_experiment_config` and `set_hydrology_parameters!`:
  use `preset`/`Config` (above). `RunSpec` has no `flux` field; the spin-up
  length is `SpinUp(years)` in the config. The never-read `log_vapor_dmc`
  switch is gone with the struct.
- `Config.modules` (never usable: `resolve` refused any value) and the unused
  `crcl` array of `CirculationWorkspace`, which has one field fewer.

### Changed

- The physics functions take the part of the configuration they read instead
  of a `PhysicsConfig`: `SWradiation!`, `LWradiation!`, `seaice!`,
  `deep_ocean!`, `advection!`, `circulation!` take a `Processes`; `hydro!`
  takes a `Processes` and a `ResolvedHydrology`; `tendencies!`, `time_loop!`,
  `qflux_correction!` and `init_model!` take a `ResolvedConfig`.
  `forcing(it, year, resolved_config)` dispatches on the scenario's parts.
  `load_cc_anomaly_jld2!` and `load_enso_anomaly_jld2!` load every anomaly
  field and take no configuration.
- `ClimateFields`: the winds split by sign are named for what they hold,
  `uclim_pos`/`uclim_neg` and `vclim_pos`/`vclim_neg` (formerly `uclim_m`,
  `uclim_p`, `vclim_m`, `vclim_p`, where `_m` held the positive part).
- `GREBClimate.derive_fields!(fields, processes)` computes every field that
  follows from the input maps (heat capacity, pressure weights, rain limit,
  deep-ocean depth, radiation-temperature offset, wind split); `init_model!`
  calls it, and it can be called again after changing a map. One land test,
  `GREBClimate.is_land(z) = z > 0`, replaces the four spellings in the
  kernels: a cell at exactly 0 m is ocean everywhere (no such cell exists in
  the dataset, so results are unchanged). The regional CO2 bands are written
  as latitudes instead of row numbers.
- Documentation: the Configuration page (formerly Physics Switches) shows
  the configuration types' own docstrings and a preset table generated from
  the presets; the API reference is split into sections.
- Source layout: the configuration lives in `src/config/`, `forcing()` in
  `src/forcing/`, `circulation.jl` in `src/physics/`.

### Changes to model results

- `:decon_2xco2` runs the MSCM physics (`mscm_hydrology()`, no moisture
  convergence) on its own 3-year spin-up, like `:decon_mean_climate`. With the
  default hydrology `humidity = :uniform` gave a non-finite run, because the
  imposed humidity is above saturation over cold, high ground and the fitted
  rain scheme rains it out at once. All processes on, the year-50 response is
  now 2.46 K (2.95 K before). Pass `hydrology = Hydrology()` and
  `processes = (moisture_convergence = true,)` for the former physics.
- `:decon_mean_climate` runs the MSCM physics (`mscm_hydrology()`, no moisture
  convergence) on the stored flux corrections. It spun up new corrections for
  every configuration before, which pulled each one back to the observed
  climate, so a switched-off process left the global mean unchanged.
- `:obliquity` and `:eccentricity` default to the table row nearest today
  (obliquity row 95, 22.5 degrees; eccentricity row 32, 0.02) instead of row 0,
  the most extreme one, whose eccentricity run failed (NaN) in year 4.
  `SolarTable(kind)` has the same default.
- Switching CO2 off (`Processes(co2 = false)`, formerly `log_co2_dmc = false`)
  sets 0 ppm in the scenario too; it applied only to the control before, and
  `:decon_mean_climate` ran its scenario at 340 ppm.
- A failed run no longer looks like a frozen planet: the 40 K floor on `Ts`
  and `Ta` turned non-finite values into exactly 40 K. They now stay NaN.
  Results of runs that stay finite are unchanged.

## [1.0.2] - 2026-09-30

### Added

- Project logo and favicon (`docs/src/assets/`), shown in the README and the
  documentation site.

### Changed

- The model's backbone (constants, config, state, tendencies, output,
  postprocess, model) now lives in `src/core/`. No API change.
- Maintainer tools are grouped by purpose; the dataset scripts moved to
  `tools/dataset/` (see `tools/README.md`).
- `tools/validation/bit_identity.jl` checks that a refactor leaves every
  output value unchanged.

### Fixed

- `greb_model!` left the climatologies replaced by the deconstruction and
  sensitivity switches (clouds, humidity, mixed layer, topography) and the
  flux corrections of one run in `fields`, so a later run on the same
  `fields` used them. The run now restores them when it returns.

## [1.0.1] - 2026-09-29

### Changed

- README and docs install the package with `Pkg.add("GREBClimate")` now that it
  is registered in the General registry.

### Changes to model results

- `RunSpec` defaults to `flux = 3`, the original GREB spin-up. With the old
  default `flux = 0` the control ran on stored flux corrections that do not fit
  the default configuration and drifted about +1.7 K over 5 years. Pass
  `flux = 0` to get the old behaviour.
- `:lanina` subtracted the ERA-Interim La Niña composite, which is already a
  cold anomaly, so it produced a warm, El Niño-like Pacific. It is now added,
  as for `:elnino` and in the Fortran (`log_exp` 241).
- `:rcp85` ran its control at 280 ppm while its scenario runs at 340 ppm, so
  the scenario carried a 60 ppm CO2 step on top of the boundary forcing. The
  control now runs at 340 ppm, as in the Fortran (`log_exp` 230).
- `:rcp85`, `:elnino` and `:lanina` added their boundary anomalies before the
  flux-correction spin-up and control run, so the control already carried the
  forced state and the returned scenario anomaly did not isolate it. They are
  now added at scenario start, as in the Fortran.
- `greb_model!` left those anomalies in `fields.Tclim`, `uclim`, `vclim`,
  `omegaclim` and `wsclim`, so a reused `fields` carried them into the next
  run and a repeated run added them twice. They are now restored when the run
  returns.

### Fixed

- `load_custom_co2_scenario` left the file open when a line was malformed.

## [1.0.0] - 2026-09-26

First release registered in the Julia General registry.

### Breaking changes

- Package renamed `GREB` -> `GREBClimate` (repository `GREBClimate.jl`); the
  General registry requires names of at least 5 characters.
- `greb_model!` refuses an unloaded `ClimateFields`: pass the value returned by
  `load_greb_jld2!` as `fields=`. Data-free runs opt in with
  `allow_uninitialized=true`.
- `create_experiment_config(:co2_double).co2_concentration` is now `340.0f0`
  (the control CO₂); `forcing` sets the scenario CO₂. Same for
  `:co2_quadruple`, `:paleo_231kyr` and `:decon_2xco2`.
- Removed: the `:a1b_enhanced` experiment (identical to `:a1b_scenario`),
  `ModelState`'s ten write-only annual-mean fields, `MonthlyAccumulator.count`
  and the unused `ε` constant.

### Changes to model results

- `:co2_double`, `:co2_quadruple`, `:paleo_231kyr` and `:decon_2xco2` showed no
  climate response: their control run used the scenario CO₂. With the Fortran's
  climatology and evaporation settings, 2×CO₂ now matches the original Fortran
  (2.758 K final-year warming after 50 years in both).
- The temperature floor `min_T_K` is 40 K instead of 233.15 K (−40 °C), which
  was clipping real polar winter cells. The control climate warms by up to
  ~0.6 K in those months; the global annual mean moves from 14.77 to 14.43 °C.
- 18 bugs found by comparison with the Fortran `greb.model.mscm.f90` were fixed.
  The largest: a stale moisture-convergence term in every temperature sub-step,
  `log_rain` having no effect, `log_eva` modes 1 and 2 duplicating mode -1, and
  `deep_ocean!` cutting heat exchange under sea ice.
- `grav` is 9.81 m/s² as in the Fortran (was 9.80665; fourth decimal).
- `:constant_topo`'s scenario runs at 680 ppm (was 550), and `:a1b_scenario`'s
  control at 280 ppm (was 298), matching the other experiments.

### Added

- **Automatic dataset download** via DataDeps.jl. `greb_data_dir()` checks an
  explicit path, `$GREB_DATA` and `greb_input_data/` first, and only then
  downloads (~353 MB, once per machine).
- **All experiments reachable from `create_experiment_config`** (42, up from
  21), with `orbital_index` and `earth_sun_distance_pct` keywords.
- **IPCC scenarios**: all four RCPs, five SSPs, a historical CO₂ hindcast
  (1850-2017) and a user-supplied CO₂ trajectory (`:custom_co2`).
- **Deconstruction experiments** `:decon_mean_climate` and `:decon_2xco2`,
  switching individual processes off.
- **Plotting toolbox** (`viz/`): maps, global-mean time series, seasonal cycle,
  Hovmöller diagram and animations, plus a simplified Pluto explorer notebook.
  `julia viz/setup.jl` sets up its environment.
- `RunSpec` for flux-correction, control and scenario lengths; explicit state
  structs (`ClimateFields`, `ModelState`, `SurfaceState`) instead of ~40
  module-level globals; keyword constructors for `ClimateFields`,
  `CirculationWorkspace` and `MonthlyAccumulator`.
- Documentation site (tutorial, input data, model overview, plots, physics
  switches, API) with doctests, and `CONTRIBUTING.md`.
- Maintainer tools: `tools/package_dataset.jl` builds the dataset archive and
  its SHA256; the `.bin` converter only converts fields the model reads.
- CI on Julia 1.10 and current, in two test shards, with code coverage,
  Aqua.jl package checks, per-kernel allocation and type-stability tests, and a
  single- vs multi-threaded bit-identity test. TagBot and CompatHelper automate
  releases and dependency bounds.

### Changed

- The README is a short landing page; the detail is in the documentation.
- The dataset shrank from 580 MB / 49 files to 439 MB / 39 files by dropping
  files no code reads. Results are unchanged.
- `forcing` is pure: the dynamic regional-CO₂ masks are built once per run.
- Source split from one 2,245-line file into topical files; tests split into
  one file per subject.

### Fixed

- The README quick start discarded the loaded data and ran on a zero
  climatology (a −40 °C world) while reporting success.
- `seaice!` returned a `Union` type; it now returns `nothing`.
- 19 of 36 exported functions had docstrings detached from their definitions.
- `tools/convert_greb_to_jld2.jl` defaulted to a non-existent `Data/input`.
- Stale references in the documentation: an old notebook name, the Julia
  version and a benchmark file that never existed.

### Performance

- About 2.3× faster per simulated year (2.7 s → 1.17 s) from running the
  temperature and humidity transport concurrently plus four `@turbo` kernels.
- `Float32` throughout: a further ~1.6×, with output within 0.01 K of the
  previous `Float64` path.
- The three flux-correction files are merged into one, ~35% faster to load.

## [0.1.0] - 2026-08-06

Initial extraction from the interactive Pluto notebook into a standard Julia
package layout (`Project.toml`, `src/`, `test/`).
