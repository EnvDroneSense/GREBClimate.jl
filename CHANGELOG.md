# Changelog

Notable changes to GREBClimate.jl, following
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Changed

- The model's backbone (constants, config, state, tendencies, output,
  postprocess, model) now lives in `src/core/`. No API change.
- Maintainer tools are grouped by purpose; the dataset scripts moved to
  `tools/dataset/` (see `tools/README.md`).

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
