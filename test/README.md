# Tests

One file per subject, `test_<subject>.jl`. `runtests.jl` includes the files
listed in the `SHARD` table of `support/selection.jl`; a file not listed there
is never run.

## Run

| Goal | Command (from the repository root) |
|:-----|:-----------------------------------|
| Tests a change needs | `julia --project=. test/support/changed.jl` (add `--run` to run them) |
| Named files | `julia --project=. -e 'using Pkg; Pkg.test(test_args=["state","io"])'` |
| A tier (`smoke` < `standard` < `full`) | `GREB_TEST_TIER=smoke julia --project=. -e 'using Pkg; Pkg.test()'` |
| Whole suite, the gate before a commit | `julia --project=. -e 'using Pkg; Pkg.test()'` |
| Whole suite in parallel processes | `julia --project=. test/support/parallel.jl [--jobs=4]` |
| Repeated runs in one session | `julia --project=. -i test/support/live.jl`, then `t("state")` |

`Pkg.test()` forces `--check-bounds=yes`; the other entry points do not.
A subset run proves only the files it ran.

## Files

| Subject | Files | Source covered |
|:--------|:------|:---------------|
| Configuration and scenarios | `test_config`, `test_processes`, `test_scenario`, `test_presets`, `test_calendar` | `src/config/`, `src/forcing/`, `src/core/calendar.jl` |
| State and bookkeeping | `test_state`, `test_output`, `test_io` | `src/core/state.jl`, `output.jl`, `postprocess.jl`, `src/io.jl`, `src/data.jl` |
| Physics kernels | `test_physics`, `test_invariants`, `test_budgets` | `src/physics/`, `src/core/tendencies.jl`, `budgets.jl` |
| Whole model | `test_model`, `test_ensemble`, `test_golden`, `test_threading` | `src/core/model.jl`, `ensemble.jl` |
| Package checks | `test_aqua`, `test_enzyme` | whole package, `ext/` |
| Maintainer tools | `test_experiment_tools`, `test_dataset_tools` (run by name only) | `tools/` |
| The test setup itself | `test_selection` | `support/selection.jl` |

The exact source-to-test map is `SOURCE_MAP` in `support/selection.jl`.

## Support

| Path | Holds |
|:-----|:------|
| `support/imports.jl` | The names the tests call besides the exports |
| `support/testutils.jl` | Fixtures: `synthetic_fields`, `run_synthetic`, `at_first_step`, file writers |
| `support/presetcapture.jl` | Captures what each preset imposes, for `test_presets` and the reference tool |
| `support/selection.jl` | `SHARD` (file, CI shard, tier, cost), `OPT_IN`, `SOURCE_MAP`, `tests_for` |
| `support/changed.jl`, `parallel.jl`, `live.jl` | The three runner scripts above |
| `support/preset_reference.jl` | Regenerates `support/data/preset_reference.jl` (on purpose only) |
| `support/data/` | Generated references: `golden_reference.jl`, `preset_reference.jl` |

## Adding a test

1. Add the testset to the file of its subject; add a file only for a new subject.
2. A new file needs a row in `SHARD` (smallest tier that includes it, rough cost in seconds).
3. A new file under `src/` or `ext/` needs a `SOURCE_MAP` row; `test_selection.jl` fails without one.
4. Tests behind `isdir(DATA_DIR)` skip in CI, which has no dataset.
