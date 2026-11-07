---
paths: ["test/**/*.jl"]
---

# Tests (layout and commands: test/README.md)

## Before running tests

- Run `julia --project=. test/support/changed.jl` (`--run` runs the files it prints). It maps `git diff` to test files through `SOURCE_MAP` in `test/support/selection.jl`; a path it does not name runs the full suite.
- The full suite is the gate before the user commits, not the loop.
- Report what ran ("ran test_state.jl, 102 passed"). A subset proves only its files.
- New test file: add a row to `SHARD` (file, CI shard, tier, cost in s). Tier: `smoke` under ~10 s; `standard` kernel, budget and preset files; `full` model runs and slow files.
- New `src/` or `ext/` file: add a `SOURCE_MAP` row (`test_selection.jl` fails without it).
- `Pkg.test()` forces `--check-bounds=yes` because the code indexes ghost-cell buffers by hand under `@inbounds`/`@turbo`. Never disable it.

## Before adding a test

Add one only if it guards public-API behaviour or a physics or budget invariant, and breaking the guarded line fails it.

- Bug fix: reproduce it through the public API (if only a hand-built internal value reaches the branch, fix the check at the public boundary and test that). Break the fixed line in a scratch copy; if an existing test fails, extend it instead.
- Grep the tests for the symbol first. The same claim in two testsets is one too many.
- A `greb_model!` run costs a simulated year: ask every question of the same config and phase of one run. Do not run the model to test what the model does not do.
- Many cases: one loop, one aggregate `@test all(...)`. Assertion count is not a metric.
- Removing a code path: delete its tests in the same commit.
- Do not test: Julia, stdlib or dependency behaviour (a JLD2 error, a `@kwdef` default); a literal against itself; a value built from the same constants in the same expression; a struct's field list or count. A cross-check between two lists that must agree is a test.
- Prefer `==` to `isapprox` where the claim is bit-identity (see the `-t 1` vs `-t 2` comparison in `test_threading.jl`).

## Fixtures (test/support/)

- `imports.jl`: names the tests call besides the exports (kernels, run loop, calendar, constants). Add a name there instead of writing `GREBClimate.name`; `_` internals keep the prefix.
- `testutils.jl`: `synthetic_fields`, `constant_fields`, `uniform_record`, `field`, `quiet`, `with_tempdir`, `with_dataset`, `gmean`, `record_times`. Input writers: `write_field`, `write_flux_corrections`, `write_boundary_anomalies`, `write_ipcc_scenarios`, `write_solar_scenarios`.
- `run_synthetic(run, config; fields, kwargs...)`: muted `greb_model!` on `synthetic_fields()` with `allow_uninitialized`, for short runs without the dataset.
- `at_first_step`: runs `greb_model!` to the first step of a phase and returns what the observer saw. Use it for what a preset decides at the start (CO2, solar table, climatology) instead of simulating a year. Not called during spin-up: pass `NoCorrections()` or `SpinUp(0)`.
- One model step: `initial_state(ini)` plus `step!(t, CO2, s, fields, r)`, which supplies the record buffers and fresh `state`, `ws`, `acc` unless given.

## Guards that need upkeep

- `test_invariants.jl` builds `all_kernels` from `names(GREBClimate; all=true)`: every `!` function in a kernel file, exported or not, needs a row in the `kernels` allocation table (equality-checked). `_` workers are skipped and covered through their wrapper; a new `_` worker with no wrapper is never checked, so give it one. The `signatures` table has no such check.
- `test_enzyme.jl` runs only where Enzyme is found (default environment or a shared `enzyme` environment); elsewhere, CI included, it skips.
