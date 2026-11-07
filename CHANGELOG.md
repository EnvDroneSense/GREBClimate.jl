# Changelog

Notable changes to GREBClimate.jl, in the spirit of
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Each version uses
these sections, in this order, and leaves out the empty ones:

| Section | Holds |
|:--|:--|
| Breaking changes | What a caller has to change: renamed or removed names, new signatures, new errors |
| Changes to model results | Every change to the numbers a run returns |
| Features | New things a user can do |
| Changes | Other changes a user can notice |
| Fixes | Wrong behaviour corrected, results unchanged |
| Performance | Faster runs or less memory |
| Documentation | The README, the docs site and the docstrings |
| Tools and tests | CI, the test suite, benchmarks and maintainer scripts |

An entry outside "Changes to model results" leaves results unchanged.

## [Unreleased]

## [2.3.0] - 2026-11-07

### Breaking changes

- The loaders and setup functions take `fields` first:
  `load_flux_corrections!(fields, dir)`,
  `load_boundary_anomaly!(fields, dir, source)`, `init_model!(fields, r)`
  (which returns `(surf, CO2_ctrl)`), `apply_co2_mask!(fields, mask)` and
  `apply_surface_mask!(fields, mask, ice_cover)`.
- `time_loop!` and `qflux_correction!` take the prognostic fields in the order
  `(Ts, Ta, To, q)`, like the physics kernels.
- The observer's view field `config` is renamed `resolved`. It holds the run's
  `ResolvedConfig`.
- The input dataset is version 2 (release `data-v2`, DataDep
  `GREB-input-data-v2`). The loaders refuse a dataset of another version with an
  `ArgumentError` that names the file and both versions. A dataset downloaded for
  2.2.x is not used and not deleted.
- The CO2 table of RCP 6.0 is `:rcp60` in the archive, and
  `load_co2_scenario(dir, :rcp6)` is an unknown-table `ArgumentError`.
- `DATA_ARCHIVE_NAME`, `DATA_SHA256` and `DATA_URL` are replaced by
  `DATA_PACKS`.

### Changes to model results

- The `ssp119`, `ssp126`, `ssp245`, `ssp460` and `ssp585` CO2 tables hold the
  historical CO2 record for 1950-2008 (dataset version 2). Before, those years
  held the RCP equivalent concentration: 310.75 against 301.92 ppm in 1950, and
  at most 15.3 ppm apart. From 2009 the tables are unchanged. Runs with an
  `ssp*` preset differ before 2009; RCP presets, `:hist` and runs without a
  table are unchanged.
- `LWradiation!`, `hydro!`, `SWradiation!` and `seaice!` no longer use
  `@turbo` (see Features). Results change at rounding level, about 3e-7
  relative.

### Features

- Optional Enzyme support: `LWradiation!`, `hydro!`, `SWradiation!` and
  `seaice!` use plain loops, so Enzyme can differentiate them. Their `log` and
  `exp` go through `fastlog` and `fastexp` (the same SLEEFPirates routines),
  whose derivative rules are in a package extension that loads with Enzyme.
  SLEEFPirates is now a direct dependency. A one-year run costs about 3 % more
  cycles.
- The loaders and `read_field` accept any `AbstractString` path.
- The dataset comes in five packs (`tar.xz`, 321 MB in all): `core` (103 MB) and
  the optional `era`, `cmip5`, `enso` and `solar`. The first `greb_data_dir()`
  asks once whether to download the whole dataset; no fetches the core pack only.
  `greb_data_dir(packs = :core)` skips the question.
- `download_greb_pack(pack; dir)` (exported) adds an optional pack, or `:all` for
  the four, to a dataset directory.
- `co2_scenario_label(dir, scenario)` (exported) returns `:co2` for the `ssp*` and
  `hist` tables and `:co2_equivalent` for the `rcp*` tables.
- Every dataset file carries `description`, `units`, `source`, `period` and
  `dataset_version`. `read_field` returns them as the new field `info`.

### Changes

- The loaders raise an `ArgumentError` for a missing file or table, an unknown
  forcing type, or an index that is not in the table.
- `read_field` raises an `ArgumentError` for a missing file (was a `jldopen`
  error). For a file of an optional pack the message names the pack.

- `viz/` has new plots: `plot_zonal`, `plot_difference` (with `to_absolute`),
  `plot_compare` (several runs as lines) and `plot_members` (one bar per ensemble
  member). `plot_timeseries` and `series` take a `mask`. New building blocks:
  `land_mask`, `band_mask`, `mask_outline`, `outline!` (`coastlines!` now uses
  them) and `zonal_mean`.
- `examples/` holds two new scripts: `quickstart.jl` (load, run, compare a
  process on and off; `--data-free` runs without the dataset) and
  `feedback_deconstruction.jl` (one `run_ensemble` that switches processes off one
  at a time and writes the warming of each to a CSV). `run_greb.jl` and
  `parameter_sweep.jl` are removed; the CO2 sweep is now
  `tools/experiments/co2_sweep.jl` with small defaults.

### Fixes

- The stand-alone `convergence!` fills the ghost cells through `to_ghosted!`,
  as the model does.

### Documentation

- The "Plots and notebook" page of the docs site is a section of four pages: setup
  and notebook, the plots, building blocks for your own plots, and how to read
  the plots. The pages show three figures, which `viz/docs_figures.jl`
  regenerates.
- A "Developer tools" section of the docs site, with a page for the dataset tools
  and pages on the benchmark harness and on the profiler that show an example
  session and how to read its output.
- The `dev` version of the docs site follows the `dev` branch instead of the
  last release. The README and CONTRIBUTING link to the `stable` version.

### Tools and tests

- `viz/notebooks/launch_pluto.jl` brings the `viz/` environment up to date on every
  start. Before, an environment set up once was never refreshed, and the notebook
  stopped loading after the package gained a dependency. The notebook and
  `viz/demo.jl` no longer copy `fields` per run.
- The preset-reference generator moved to `test/support/preset_reference.jl` (was
  `tools/validation/`) and runs again; it failed to load. `bit_identity.jl` and
  `budget.jl` can be `include`d without exiting the session.
- The dataset tools read the field metadata from `tools/dataset/fields.toml`.
  `convert_greb_to_jld2.jl` stops when a field fails, prints a summary per
  section and no longer reads `--all` as the input directory;
  `package_dataset.jl` builds the five packs, requires every combined file and
  checks the metadata.
- `test/test_dataset_tools.jl` tests these tools and that their file list matches
  the loaders. It runs only by name (`Pkg.test(test_args=["dataset_tools"])`),
  not in the full suite.
- CI also runs on pushes to `dev`, and a CI job checks the formatting against
  `.JuliaFormatter.toml`.
- The benchmark harness times with Chairmarks. Its new `compare` mode runs the
  variants in a shuffled order each round and checks itself against a second
  copy of the baseline.
- `benchmark/run_benchmarks.jl` times each run once (it ran the model twice per
  run and copied the climatology each time), and `year` is the `years` mode with
  one control year. The `alloc` mode is removed: `test/test_invariants.jl` checks
  the same budget. A mistyped flag, a reps argument a mode does not take, an
  unknown mode and a missing dataset now stop the run with exit code 1.
  `threads` reports the median, as the other modes do. CI also loads the harness.
- `benchmark/profile.jl` writes one `report.md` per run (Markdown tables that
  also read as plain text) in place of `header.txt`, `category.txt`, `owned.txt`,
  `allocs.txt` and `dispatch.txt`, and prints the top rows of the same tables.
  `step` and `setup` also write `counts.tsv`, which `compare` reads instead of
  parsing text reports; rows below a run's cut-off no longer rank first, so the
  `?` marker is gone, and runs made before this change must be profiled again.
  `tree.txt` and `flat.txt` now start at `greb_model!` and not at the profiler's
  own frames, the script returns an exit code and rejects extra arguments, and
  the helpers it shares with `run_benchmarks.jl` moved to `common.jl`.
- The Fortran comparison (`benchmark/fortran/run_fortran_comparison.jl`, with
  `peak_memory.ps1`) is now in the repository. The Julia side is a function of
  the script, started in a child process, in place of a generated text script, so
  a renamed package function fails when the script loads. A failed gate,
  missing tooling or a failed self-check exits with code 1. `verify` also fails
  when the Fortran anomaly is below 0.1 K, the input link is recreated when
  `GREB_FORTRAN_DIR` changes, and the Windows-only linker flag is no longer
  passed on Linux. `compare` now prints the machine state like the other
  harnesses and appends its rows to `benchmark/results/` (`record_result` takes
  a `threads` keyword for this). The script uses the helpers in `common.jl`, and
  the developer docs have a page "Fortran comparison".
- CI checks the formatting of `tools`, `ext`, `examples` and `viz` as well as
  `src`, `test` and `benchmark`, runs the `viz` tests, loads the profiler, and
  runs the Fortran harness self-test. Those folders were formatted.
- The developer docs have pages for the diagnostics, experiments and validation
  tools, for testing and for the code conventions, and the tutorial has an
  Examples page with the output of both example scripts. `CONTRIBUTING.md` is
  shorter and links to them; the dataset size is given as 439 MB unpacked and a
  321 MB download everywhere.
- The test that every physics kernel has an allocation budget checked nothing
  since 2.1.0, because it listed only exported names. It now finds unexported
  kernels, fails when it finds none, and covers `to_ghosted!` and
  `refresh_ghosts!`.

## [2.2.0] - 2026-10-09

### Breaking changes

- One time value, `ModelTime(phase, start_year, step)`, says where a run is;
  the day, month and year are derived from it:
  - `greb_model!` returns `(ctrl, scnr, ctrl_time, scnr_time, scnr_anomaly)`.
    `ctrl_time` and `scnr_time` hold the `(year, month)` of each record, and
    `scnr_anomaly` says whether `scnr` is an anomaly against the control.
  - `forcing(t::ModelTime, r)` replaces `forcing(it, year, r)`.
  - `scenario_anomalies(scnr, scnr_time, ctrl_clim)` takes the record times and
    no longer assumes that the first record is January.
  - The observer's view has `time`, `phase`, `step`, `year` and `step_of_year`
    in place of `phase`, `it`, `year` and `ityr`. `phase` is
    `GREBClimate.control` or `GREBClimate.scenario` (was `:ctrl` or `:scnr`),
    and `RangeCheck`'s `first` has `step` in place of `it`.
  - `time_loop!`, `output!` and `diagnostics!` take the time value and return
    `nothing`.
- `TimeState` is removed. The physics kernels and `tendencies!` take the index
  of the climatology slice, an `Int`; `time_loop!` and `qflux_correction!` lose
  the argument.
- `MonthlyAccumulator`'s fields take the names of `MonthlyRecord` (`Tmm` is
  `Ts`, `apmm` is `albedo`, `olrmm` is `olr`, and so on), declared once in
  `GREBClimate.OUTPUT_FIELDS`. Its keyword constructor is gone.
- `ClimateFields` no longer has the ten anomaly arrays (`Ts_clim_anom_enso`,
  `Ts_clim_anom_cc` and the others). A `BoundaryAnomaly` source is loaded into
  `fields.boundary_anomaly`, which is `nothing` until then, and
  `anom_cc_source` and `anom_enso_source` become `boundary_anomaly_source`.
- `notebooks/` moved to `viz/notebooks/`: launch with
  `julia viz/notebooks/launch_pluto.jl`.

### Changes to model results

- `:a1b` (`A1BRamp`) holds 700 ppm after 2100, as the Fortran does. It fell
  back to 340 ppm in 2101. Runs that end by 2100 are unchanged.

### Features

- The plotting toolbox (`viz/`) reads the year and month of each record from
  `ctrl_time` and `scnr_time`, and whether the scenario is an anomaly from
  `scnr_anomaly` (it used to guess from the temperature). Animation frames are
  titled `Jul 1950`. `annual`, `seasonal_cycle` and `map_frames` take the times
  as an optional argument.

### Changes

- Calendar constants are written out: `ndt_days`, `ndays_yr`, `cjday_mon` and
  `jday_mon_cumsum` are `steps_per_day`, `days_per_year`, `days_in_month` and
  `last_day_of_month`.

### Fixes

- `run_ensemble` kept the copies of the fields alive after it returned, 2 to
  3 GiB more than needed.

### Performance

- A `ClimateFields` is 135 MB smaller, also per ensemble member (see the
  anomaly arrays above).
- A `Stored` run on a `fields` that already holds the corrections from the same
  file no longer saves, reloads and restores them: 4 to 7 % of a one-year run.
  `ClimateFields.flux_source` names that file; set it to `""` after writing into
  the correction arrays.

### Tools and tests

- `test/parallel.jl` runs the test files in several processes: about 90 s
  instead of about 295 s on a 14-thread machine.
- `benchmark/profile.jl step` samples 25 years by default (was 50).
- `examples/parameter_sweep.jl` runs its CO2 levels side by side with
  `run_ensemble`. `tools/experiments/` holds helpers for experiment scripts and
  `sensitivity.jl`, which prints the effective climate sensitivity of an
  abrupt-CO2 run.

## [2.1.0] - 2026-10-05

### Breaking changes

- The loaders no longer carry the file format in their names, and the two
  anomaly loaders are one function:

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

- `ClimateFields` fields separate their words:

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
- `CirculationWorkspace` is `ModelWorkspace`, and its 24 fields lose the `_buf`
  and `_out` suffixes (`ws.Q_lat_buf` is `ws.Q_lat`).
- Other renames: the `Processes` switches `vapour_diffusion` and
  `vapour_advection` are `vapor_diffusion` and `vapor_advection`;
  `build_monthly_climatology` is `monthly_climatology`,
  `apply_scenario_anomalies` is `scenario_anomalies`,
  `compute_annual_ice_climatology` is `ice_climatology` and
  `apply_dynamic_co2_mask!` is `apply_surface_mask!`.
- Scenario fields say what they hold: `SSTOffset.offset` (was `K`),
  `SolarConstant.offset` (`dW`), `EarthSunDistance.percent` (`pct`, also the
  preset keyword), `SeasonalCO2.in_season` and `.out_of_season` (`inside`,
  `outside`). Positional construction is unchanged.
- The kernels and loop functions are no longer exported: `SWradiation!`,
  `LWradiation!`, `hydro!`, `convergence!`, `seaice!`, `deep_ocean!`,
  `diffusion!`, `advection!`, `circulation!`, `tendencies!`, `time_loop!`,
  `output!`, `diagnostics!`, `qflux_correction!` and `init_model!`. Call them
  as `GREBClimate.name`.
- Arguments that were never read are gone: `diagnostics!(year, surf, state,
  timestate)` and `diffusion!(T1, h_scl, fields, ws)`.
- `ModelWorkspace` and `MonthlyAccumulator` are immutable; their arrays are
  still written in place.
- The RCP6.0 CO2 table is `CO2Table(:rcp60)`, the name of its preset;
  `:rcp6` raises an `ArgumentError`.
- `load_climatology` and `load_flux_corrections!` raise an `ArgumentError` when
  `climatology/flux_corrections.jld2` or one of its tables is missing, instead
  of filling zeros with a warning. Pass `corrections = false` for a dataset
  without the file.

### Features

- `run_ensemble([reduce,] run, configs; fields, jld2_dir, ntasks, logger)` runs
  many configurations side by side, each giving exactly the result of running
  it alone.
- `global_mean(field)`: the area-weighted global mean. The examples use it;
  their plain means read 8 to 11 K too cold.
- `GREBClimate.RangeCheck()`: an observer that records the first step at which
  `Ts`, `Ta`, `To` or `q` leaves a physical range.
- Calendar functions `step_of_year`, `day_of_year`, `month_of_step` and
  `decimal_year`.
- `load_climatology(dir; corrections = false)` loads a dataset without the
  stored flux corrections.

### Changes

- A run and the loaders report through the logger (`@info`) instead of
  `println`: silence them with `with_logger(NullLogger())`.
- The yearly line says what it shows
  (`1970: Ts global mean 13.82 °C; 178 E 9 N 26.79; 58 E 51 N 4.85`), a spin-up
  year reads `spin-up year 1: ...`, and a run without a scenario no longer
  reports a scenario of 0 years.
- Internal constants in the model's formulas are lower case (`is_polar`,
  `polar_diff_time2`, `ΔT_air_factor`), and some are named for what they are:
  `pi_f32`, `convergence_factor`, `q_to_mm_per_day`, `solar_percent` and
  `emissivity_fit`. No exported name is affected.

### Performance

- A `BoundaryAnomaly` scenario (`:rcp85_boundary`, `:elnino`, `:lanina`) reads
  its files once per `ClimateFields` and source.
- The calendar is computed instead of looked up: 2.3 MB less memory.
- The precompile workload runs a scenario year too, so the first scenario run
  of a session compiles less.

## [2.0.1] - 2026-10-04

### Changes to model results

- `:regional_co2_nh`, `:regional_co2_sh`, `:regional_co2_tropics` and
  `:regional_co2_extratropics` (any `LatitudeMask`): the mask was already
  active in the spin-up and the control, which ran on half the control CO2
  outside the band. It now applies from the start of the scenario, as in the
  Fortran. The response outside the band is much smaller than before.

### Changes

- The `observer` keyword of `greb_model!` and `GREBClimate.BudgetCheck` are no
  longer experimental; the observer's `view` is part of the interface.
- `CirculationWorkspace` loses seven fields that nothing read: `qs`, `rq`,
  `Tskin`, `ws_base`, `cE_buf`, `temp_buf` and `a_atmos_buf`.

### Fixes

- `corrections = Stored()` without `jld2_dir` ran without corrections. It now
  uses the corrections already in `fields`; with `jld2_dir`, a missing
  corrections file is an `ArgumentError`.
- `load_greb_jld2!` with an unknown `dataset` silently loaded NCEP; it now
  raises an `ArgumentError`.
- A `SurfaceMask` scenario with `RunSpec(ctrl = 0)` built its mask from an ice
  cover of zero; it now raises an `ArgumentError`.
- Runs in parallel tasks could fail with a JLD2 error when they read dataset
  files at the same time; the loaders now lock each file open.
- An orbital `index` that is not in the table, or a solar table of the wrong
  shape, raise an error naming the problem instead of an `@assert`.

### Performance

- An `SSTOffset` scenario (`:sst_plus1`) no longer copies the temperature
  climatology at every step.

### Documentation

- Docstrings corrected for `CO2Table`, `Hydrology` (`:skin_gust`),
  `ClimateFields`, `tendencies!` and `mscm_hydrology`.

## [2.0.0] - 2026-10-01

### Breaking changes

- A new configuration API replaces `PhysicsConfig`, `create_experiment_config`
  and `set_hydrology_parameters!` (see Features). `RunSpec` has no `flux`
  field: the spin-up length is `SpinUp(years)` in the config.
- Preset renames: `:rcp85` is `:rcp85_boundary`, `:co2_step` is
  `:co2_abrupt_reverse` and `:a1b_scenario` is `:a1b`. `:constant_topo` has no
  preset; use `preset(:co2_double; processes = (topography = :flat,),
  corrections = Stored())`.
- The physics functions take the part of the configuration they read: a
  `Processes`, a `ResolvedHydrology` (`hydro!`) or a `ResolvedConfig`
  (`tendencies!`, `time_loop!`, `qflux_correction!`, `init_model!`).
  `load_cc_anomaly_jld2!` and `load_enso_anomaly_jld2!` take no configuration.
- `MonthlyRecord` has two more fields, `olr` and `lwdown`; code that builds one
  by hand must add them.
- The split winds are named for what they hold: `uclim_pos`, `uclim_neg`,
  `vclim_pos`, `vclim_neg` (were `uclim_m`, `uclim_p`, `vclim_m`, `vclim_p`).
- Removed: `Config.modules` (never usable), the `log_vapor_dmc` switch (never
  read) and the unused `crcl` array of `CirculationWorkspace`.

### Changes to model results

- `:decon_2xco2` runs the MSCM physics (`mscm_hydrology()`, no moisture
  convergence) on its own 3-year spin-up. With the default hydrology,
  `humidity = :uniform` gave a non-finite run. All processes on, the year-50
  response is 2.46 K (was 2.95 K).
- `:decon_mean_climate` runs the MSCM physics on the stored flux corrections.
  It used to spin up for every configuration, which pulled the climate back to
  the observed one, so a switched-off process left the global mean unchanged.
- `:obliquity` and `:eccentricity` default to the table row nearest today
  (obliquity row 95, eccentricity row 32) instead of row 0, whose eccentricity
  run failed in year 4. `SolarTable(kind)` has the same default.
- Switching CO2 off (`Processes(co2 = false)`) sets 0 ppm in the scenario too,
  not only in the control.
- A failed run no longer looks like a frozen planet: non-finite `Ts` and `Ta`
  stay NaN instead of being floored to 40 K.

### Features

- The configuration API: a `Config` holds a `Scenario` (CO2 path and mask,
  sunlight, surface forcing, control CO2, start year, output form), a
  `Processes` set (one option per process, each with a stated meaning when
  off), a `Hydrology` scheme (`mscm_hydrology()` for the original) and
  `Corrections` (`SpinUp(years)`, `Stored()`, `NoCorrections()`).
  `preset(name)` builds every experiment with the same results as before, and
  `preset(name; processes = (ocean = :none,))` changes single options. Scenario
  parts combine freely. `resolve(config)` returns a `ResolvedConfig` with the
  tables and coefficients the configuration refers to.
- `:rcp85`: RCP8.5 CO2 from the dataset's table, like the other IPCC presets.
- `MonthlyRecord` fields `olr` (longwave leaving to space, positive upward) and
  `lwdown` (longwave reaching the surface, positive downward), in W/m2.
- `greb_model!(...; observer = f)` (experimental): `f(point, view)` is called
  before and after every step with the model's state and flows.
  `GREBClimate.BudgetCheck()` is an observer that checks each store changes by
  the sum of its flows.
- `GREBClimate.derive_fields!(fields, processes)` recomputes every field that
  follows from the input maps, for use after changing a map.

### Changes

- One land test, `GREBClimate.is_land(z) = z > 0`, replaces four spellings: a
  cell at exactly 0 m is ocean everywhere (the dataset has none).
- The regional CO2 bands are written as latitudes instead of row numbers.
- The configuration lives in `src/config/`, `forcing()` in `src/forcing/` and
  `circulation.jl` in `src/physics/`.

### Documentation

- The Configuration page (formerly Physics Switches) shows the configuration
  types' docstrings and a preset table generated from the presets. The API
  reference is split into sections.

### Tools and tests

- A regression test that the MSCM configuration reproduces the MSCM 2xCO2
  response (year 1 global mean 0.5946 K).
- `tools/diagnostics/budget.jl` prints the result of a `BudgetCheck` and the
  global-mean flows of a run.

## [1.0.2] - 2026-09-30

### Changes

- The model's backbone lives in `src/core/`. No API change.

### Fixes

- `greb_model!` left the climatologies changed by the deconstruction and
  sensitivity switches, and one run's flux corrections, in `fields`, so a later
  run on the same `fields` used them. The run now restores them.

### Documentation

- Project logo and favicon in the README and on the docs site.

### Tools and tests

- Maintainer tools are grouped by purpose (see `tools/README.md`).
- `tools/validation/bit_identity.jl` checks that a refactor leaves every output
  value unchanged.

## [1.0.1] - 2026-09-29

### Changes to model results

- `RunSpec` defaults to `flux = 3`, the original GREB spin-up. With `flux = 0`
  the control drifted about +1.7 K over 5 years on stored corrections that do
  not fit the default configuration.
- `:lanina` subtracted the La Niña composite, which is already a cold anomaly,
  and produced an El Niño-like Pacific. It is now added, as in the Fortran.
- `:rcp85` ran its control at 280 ppm and its scenario at 340 ppm; the control
  now runs at 340 ppm, as in the Fortran.
- `:rcp85`, `:elnino` and `:lanina` added their anomalies before the spin-up
  and control, so the returned anomaly did not isolate the forcing. They are
  now added at scenario start, as in the Fortran.
- `greb_model!` left those anomalies in `fields`, so a reused `fields` added
  them twice. They are now removed when the run returns.

### Fixes

- `load_custom_co2_scenario` left the file open when a line was malformed.

### Documentation

- README and docs install with `Pkg.add("GREBClimate")`.

## [1.0.0] - 2026-09-26

First release registered in the Julia General registry.

### Breaking changes

- Package renamed `GREB` -> `GREBClimate` (repository `GREBClimate.jl`).
- `greb_model!` refuses an unloaded `ClimateFields`; data-free runs opt in with
  `allow_uninitialized=true`.
- `create_experiment_config(:co2_double).co2_concentration` is the control CO2
  (`340.0f0`), and `forcing` sets the scenario CO2. Same for `:co2_quadruple`,
  `:paleo_231kyr` and `:decon_2xco2`.
- Removed: the `:a1b_enhanced` experiment (identical to `:a1b_scenario`),
  ten write-only fields of `ModelState`, `MonthlyAccumulator.count` and the
  unused `ε` constant.

### Changes to model results

- `:co2_double`, `:co2_quadruple`, `:paleo_231kyr` and `:decon_2xco2` showed no
  response, because their control used the scenario CO2. 2xCO2 now matches the
  Fortran (2.758 K after 50 years in both).
- The temperature floor `min_T_K` is 40 K instead of 233.15 K, which clipped
  real polar winter cells. The global annual mean moves from 14.77 to 14.43 °C.
- 18 bugs found by comparison with the Fortran were fixed. The largest: a stale
  moisture-convergence term in every temperature sub-step, `log_rain` having no
  effect, two evaporation modes duplicating a third, and `deep_ocean!` cutting
  heat exchange under sea ice.
- `grav` is 9.81 m/s2 as in the Fortran (was 9.80665).
- `:constant_topo`'s scenario runs at 680 ppm (was 550), and `:a1b_scenario`'s
  control at 280 ppm (was 298).

### Features

- Automatic dataset download via DataDeps.jl, after an explicit path,
  `$GREB_DATA` and `greb_input_data/` are checked.
- All 42 experiments are reachable from `create_experiment_config` (was 21):
  all four RCPs, five SSPs, a historical CO2 hindcast (1850-2017), a
  user-supplied CO2 path (`:custom_co2`), orbital experiments and the
  deconstruction experiments `:decon_mean_climate` and `:decon_2xco2`.
- A plotting toolbox (`viz/`) with maps, time series, seasonal cycle,
  Hovmöller diagram and animations, and a Pluto explorer notebook.
- `RunSpec`, and explicit state structs (`ClimateFields`, `ModelState`,
  `SurfaceState`) instead of about 40 module-level globals.

### Changes

- The dataset dropped files no code reads: 439 MB and 39 files, from 580 MB and
  49.
- `forcing` is pure: the regional CO2 masks are built once per run.
- The source is split into topical files and the tests into one file per
  subject.

### Fixes

- The README quick start ran on a zero climatology while reporting success.
- `seaice!` returned a `Union` type; it now returns `nothing`.
- `tools/convert_greb_to_jld2.jl` defaulted to a missing `Data/input`.

### Performance

- About 2.3x faster per simulated year (2.7 s to 1.17 s) from running the
  temperature and humidity transport concurrently and four `@turbo` kernels.
- `Float32` throughout: a further 1.6x, with output within 0.01 K of the
  `Float64` path.
- The three flux-correction files are merged into one, about 35 % faster to
  load.

### Documentation

- A documentation site (tutorial, input data, model overview, plots, physics
  switches, API) with doctests, a short README and `CONTRIBUTING.md`.
- 19 of 36 exported functions had docstrings detached from their definitions,
  and the docs had stale references.

### Tools and tests

- CI on Julia 1.10 and current in two test shards, with coverage, Aqua.jl,
  per-kernel allocation and type-stability tests and a threaded bit-identity
  test. TagBot and CompatHelper automate releases and dependency bounds.
- `tools/package_dataset.jl` builds the dataset archive and its SHA256; the
  `.bin` converter converts only fields the model reads.

## [0.1.0] - 2026-08-06

Initial extraction from the interactive Pluto notebook into a standard Julia
package layout (`Project.toml`, `src/`, `test/`).
