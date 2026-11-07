# Testing

The test suite is the correctness gate. `test/README.md` has the file map and the support files; this page says how to run the suite and what a new test has to do.

| Goal | Command (from the repository root) |
|:-----|:-----------------------------------|
| Whole suite | `julia --project=. -e 'using Pkg; Pkg.test()'` |
| Whole suite, in parallel processes | `julia --project=. test/support/parallel.jl [--jobs=4]` |
| One CI shard | `GREB_TEST_SHARD=light julia --project=. -e 'using Pkg; Pkg.test()'` (or `heavy`) |
| Named files | `julia --project=. -e 'using Pkg; Pkg.test(test_args=["state","io"])'` |
| A tier | `GREB_TEST_TIER=smoke julia --project=. -e 'using Pkg; Pkg.test()'` (`smoke`, `standard`, `full`; each includes the ones before) |
| The files a change needs | `julia --project=. test/support/changed.jl` (add `--run` to run them) |
| Repeated runs in one session | `julia --project=. -i test/support/live.jl`, then `t("state")` or `t(:smoke)` |

The suite has two shards, grouped by measured run time and only roughly balanced. CI runs them as separate jobs. The parallel runner takes about 100 s instead of about 295 s and needs roughly 1 GB of memory per process. Each process runs `Pkg.test`, so bounds checking stays on.

A run of a subset proves only the files it ran. Run the whole suite before you open a pull request.

## Writing a test

| Rule | Why |
|:-----|:----|
| Add a new test file to the `SHARD` table in `test/support/selection.jl` | Nothing globs the directory, so a file that is not listed never runs |
| Pass `allow_download=false` to every `greb_data_dir` call in a test | Without it the test passes locally, where a dataset exists, and fails in CI |
| Do not run the model to test something the model does not do | Most of the suite's cost is simulated years. Call `init_model!`, a loader or a single kernel directly |
| Use `at_first_step()` for what a run decides at its start (the CO2 or solar table, the climatology) | It runs `greb_model!` to the first step, returns what the observer saw and stops. Pass `NoCorrections()`: the observer is not called during the spin-up |
| Break the line a test claims to pin and check that the test fails | A test that cannot fail guards nothing. When you delete a test, check that another one fails too |
| Use fixtures from `test/support/testutils.jl` | `quiet()`, `with_tempdir()`, `synthetic_fields()`, `constant_fields()`, `uniform_record()`, `at_first_step()`, `DATA_DIR` and the grid constants |

## What CI cannot check

CI has no dataset, so every test behind `isdir(DATA_DIR)` skips there: the golden regression, the MSCM, deconstruction and orbital-table tests, and the threaded-against-serial comparison on the real dataset. A green CI says nothing about them. Run the whole suite locally with the dataset before you open a pull request.

| File | Rule |
|:-----|:-----|
| `test/test_golden.jl` | A regression lock on model output. If your change moves those numbers, it is either a bug or a deliberate physics change that the pull request must say |
| `test/test_presets.jl` | Compares every preset with `test/support/data/preset_reference.jl`. After a deliberate preset change, regenerate the file in the same commit with `julia --project=. test/support/preset_reference.jl` and say so in the pull request. Never edit it by hand |
| `test/test_invariants.jl` | Checks the allocation budget and the return type of every physics kernel. A new kernel needs a row in both tables. Only the allocation table is checked for completeness, so add the return-type row yourself |
