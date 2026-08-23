# Changelog

Notable changes to GREBClimate.jl, in roughly chronological order. Loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) conventions;
since the package hasn't had its first registered release yet, entries are
grouped by development phase rather than version number. Each entry below records what
changed; the pass-by-pass discovery narrative behind them is kept in the
maintainers' working notes rather than in this file.

## [Unreleased] - toward v1.0.0 / first registration

### Renamed
- Package `GREB` → `GREBClimate`; repository `GREB-julia` → `GREBClimate.jl`
  (the General registry's AutoMerge requires package names of at least 5
  characters).

### Added
- **Automatic dataset download.** The JLD2 input dataset is now fetched and
  cached on first use via [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl).
  The new exported `greb_data_dir()` resolves an explicit path, then
  `$GREB_DATA`, then a local `greb_input_data/`, and only then downloads - so
  existing setups keep working untouched and offline use is unaffected. Pass
  `allow_download=false` to stop before the network; the test suite and the
  benchmarks do exactly that, so neither can pull 353 MB as a side effect.
  Set `DATADEPS_ALWAYS_ACCEPT=true` for non-interactive sessions.
- `tools/package_dataset.jl` builds the published dataset archive
  reproducibly and prints its SHA256, validating the tree against the
  converter's allowlist first so a stray field cannot enlarge the download.
- `greb_model!` now refuses an uninitialized `ClimateFields` (all-zero
  climatology) rather than silently simulating a physically meaningless
  ~-40 °C world. `ClimateFields` carries a `loaded` flag set by
  `load_greb_jld2!`; genuinely data-free runs - package precompilation and
  config/scenario-plumbing tests - opt in with `allow_uninitialized=true`.
- `tools/convert_greb_to_jld2.jl` filters on an explicit `MODEL_FIELD_NAMES`
  allowlist, so it no longer emits 11 `.jld2` files (~148 MB) that nothing
  reads. `--all` restores the previous convert-everything behaviour. A test
  asserts the allowlist and `src/io.jl`'s loads agree in both directions.
- Documenter.jl site - tutorial, API reference, and a physics-switches guide
  - deployed to GitHub Pages.
- CI: matrix over Julia `1.10`/`1`, split into `light`/`heavy` test shards.
- A benchmark suite (`benchmark/run_benchmarks.jl`) with `year`/`stages`/
  `threads`/`alloc` modes.
- `ClimateFields`/`ModelState`/`SurfaceState` structs, replacing ~40 mutable
  module-level globals with explicit, passed-in state.
- `RunSpec` for flux-correction/control/scenario run lengths.
- All four IPCC RCP scenarios and five SSP scenarios, a historical CO₂
  (1850–2017) hindcast, and a user-supplied CO₂ trajectory (`:custom_co2`).
- "Deconstruct" experiment presets (`:decon_mean_climate`, `:decon_2xco2`)
  toggling individual feedback processes.
- 3-way threaded physics: the temperature and humidity `circulation!` calls
  run concurrently via `Threads.@spawn`.
- `Float32` compute path throughout, matching the JLD2 input data's native
  precision.
- **`create_experiment_config` now reaches every experiment the model
  dispatches on** - 42, up from 21. The 22 that previously required
  hand-building a `PhysicsConfig` (all eight `regional_co2_*`, `:obliquity`,
  `:eccentricity`, `:earth_sun_distance`, `:co2_10x`, `:co2_half`, `:co2_zero`,
  `:co2_sine_wave`, `:co2_step`, `:solar_cycle_11yr`, `:a1b_scenario`,
  `:sst_plus1`, and the two crossed paleo/modern presets) are now first-class.
  New `orbital_index` and `earth_sun_distance_pct` keywords plumb the
  parameters those experiments need. An unknown symbol errors with the valid
  list instead of a bare message.
- Keyword constructors for `ClimateFields`, `CirculationWorkspace` and
  `MonthlyAccumulator` (`Base.@kwdef`), e.g. `ClimateFields(loaded = true)`.
  The positional and no-argument forms are unchanged. `ClimateFields`'s 39
  fields were previously matched to 39 bare `zeros(...)` calls by ordinal
  position alone, where a reordering would have silently assigned the wrong
  array to the wrong field.

### Removed
- **The `:a1b_enhanced` experiment.** Its CO₂ ramp was byte-identical to
  `:a1b_scenario`'s; use that instead.
- `ModelState`'s ten annual-mean accumulator fields (`Tamn`, `Tomn`, `qmn`,
  `amn`, `swmn`, `lwmn`, `qlatmn`, `qsensmn`, `ftmn`, `fqmn`) and
  `MonthlyAccumulator.count`. All were write-only: accumulated, divided and
  zeroed inside a single `diagnostics!` call without ever being read, so they
  could not be used as an output path. Model output is the
  `Vector{MonthlyRecord}` `greb_model!` returns. `ModelState` keeps `Tsmn`,
  which backs the printed annual progress line.
- The unused `ε` (IR emissivity) constant.

### Fixed
- **`min_T_K` was clamping legitimate polar temperatures.** The floor was
  233.15 K (−40 °C), cold enough to silently truncate real Antarctic and
  Siberian winter cells; it is now 40 K, a pure numerical-stability floor.
  **This changes results**: the control climate warms by up to ~0.6 K in the
  affected months, and the printed global annual mean moves from 14.77 °C to
  14.43 °C. The golden-regression snapshot was regenerated to match. A
  data-free run now reports a 40 K world rather than a 233.15 K one.
- `grav` corrected from 9.80665 to 9.81 m/s², matching the reference Fortran
  GREB. Affects results at the fourth decimal.
- `:constant_topo`'s scenario CO₂ was 550 ppm, an arbitrary value; it is now
  680 ppm (2×340), consistent with the other doubling experiments.
- `:a1b_scenario`'s control baseline moved from 298 ppm to 280 ppm, matching
  every other IPCC-style scenario.
- `forcing` is now pure. The `:regional_co2_ocean`/`:regional_co2_land_ice`
  masks were built inside it behind an `it == 1` guard, making a per-timestep
  function statefully non-idempotent; they are now computed once per run
  between the control and scenario phases. Same masks, same ordering - the
  control run still sees an all-ones mask.
- **The multithreaded circulation path was never tested.** `tendencies!` runs
  `circulation!(Ta)`/`circulation!(q)` concurrently only when
  `Threads.nthreads() > 1`, and the test process was single-threaded, so that
  branch had only ever been validated by benchmarking. A heavy-shard test now
  spawns `-t 1` and `-t 2` subprocesses and asserts bit-identical monthly
  means, and CI sets `JULIA_NUM_THREADS=2`.
- **README quick-start produced wrong results.** The documented first run
  called `load_greb_jld2!` and discarded its return value, so the model ran on
  a zero climatology and reported a global-mean surface temperature of
  233 K (-40 °C) instead of 277 K (14.8 °C) - while printing
  `✅ All GREB data loaded successfully` and passing `all(isfinite, Ts)`. The
  snippet now threads `fields` through, and the API refuses the mistake.
- A second instance of the same pattern in README §Loading Data.
- `tools/convert_greb_to_jld2.jl` defaulted its input directory to
  `Data/input`, which does not exist - the layout is flat `Data/` plus
  `Data/solar_forcing_scenarios/` - so the documented no-argument invocation
  always failed. Corrected, with an actionable error when the directory is
  missing.
- `docs/src/index.md` claimed Julia 1.9; `Project.toml` requires 1.10.
- `docs/src/tutorial.md` referenced a `claude/BENCHMARKS.md` that never existed.
18 correctness bugs found by direct comparison against the Fortran
reference (`greb.model.mscm.f90`); the ones with the widest-reaching impact:

- `circulation!` left a stale moisture-convergence term in every
  temperature sub-step of every default run.
- `set_hydrology_parameters!` wrote to disconnected module globals instead
  of the config struct - the `log_rain` switch had zero effect on any run.
- `hydro!`'s `log_eva` modes `1`/`2` silently duplicated mode `-1`'s formula
  instead of using their own coefficients.
- `deep_ocean!` cut off ocean-atmosphere heat exchange under sea ice and in
  high-latitude winters.
- 19 of 36 exported functions had docstrings silently detached from their
  definitions by a comment-placement pattern that breaks Julia's `@doc`
  binding - found while wiring up the Documenter.jl site.

The remaining 13 fixes, plus one investigated-and-reverted finding, are
recorded in the maintainers' working notes.

### Performance
- ~2.31× faster per simulated year (2.7s → 1.17s) from 3-way threading of
  `tendencies!` combined with 4 `@turbo`-rewritten hot-loop kernels.
- `Float32` throughout: a further ~1.6× on top of threading, with output
  validated against the previous `Float64` path to well under 0.01 K.
- Flux-correction files (`Tsurf`/`vapour`/`Tocean` corrections) merged into
  a single `flux_corrections.jld2` - ~35% faster to load, no size penalty.

### Changed
- The dataset shrank from 580 MB / 49 files to **439 MB / 39 files**: 11 files
  that no code reads were removed, mostly CMIP5 `.new` variants of fields the
  model reads in their non-`.new` form. The official MSCM Fortran GREB opens
  the non-`.new` names, so this changes no results.
- User-facing data documentation now describes the automatic download and the
  paths that bypass it. `DATA_README.md` and the `.bin` converter are labelled
  as maintainer tooling, which is what they are - the raw inputs are collated
  from several upstream sources and are not redistributed.
- `src/GREB.jl`, originally a single 2,245-line file, split into topical
  files (`constants`, `config`, `state`, `io`, `physics/`, `circulation`,
  `tendencies`, `output`, `postprocess`, `model`).
- Test suite split into `light`/`heavy` shards (`GREB_TEST_SHARD` env var),
  matching the CI matrix.
- `test/runtests.jl`, a single 1,407-line file, split into one file per subject
  (`test_config`, `test_state`, `test_output`, `test_physics`, `test_io`,
  `test_threading`, `test_model`, `test_golden`) plus shared fixtures in
  `test/testutils.jl`. `runtests.jl` is now a 30-line table mapping each file to
  its shard, so the split is declarative rather than two hand-maintained
  function bodies.
- The suite runs in **1m56s, down from 3m46s**, with two more assertions than
  before (1193 vs 1191). The saving came from not running the model to test
  things the model does not do: the `co2_part` reset lives in `init_model!`, the
  scenario tables are read by `load_co2_scenario_jld2`/`load_custom_co2_scenario`/
  `load_solar_forcing_jld2`, and three of the four `log_eva` branches are
  reachable by calling `hydro!` directly - 22 simulated years became 8. The
  threaded-vs-serial subprocess test now uses synthetic fields instead of
  loading the 390 MB dataset twice (50s → 9.5s), which also means it runs in CI
  rather than skipping. Shards are now balanced by measured runtime (54s/72s,
  previously 17s/176s), cutting CI wall clock from 176s to 72s.

## [0.1.0] - 2026-08-06
Initial extraction from the interactive Pluto notebook into a standard Julia
package layout (`Project.toml`, `src/`, `test/`).
