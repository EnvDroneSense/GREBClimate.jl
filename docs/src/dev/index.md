# Developer tools

The scripts in `tools/` and the timing harness in `benchmark/` are for maintainers. Nobody needs them to use the package. Run each from the repository root, the tools with `julia --project=. tools/<folder>/<script>.jl`. This section has one page per group of tools; every group has a page here, and `tools/README.md` keeps one row per script.

| Folder | Purpose | Page |
|:-------|:--------|:-----|
| `tools/dataset/` | Build and publish the input dataset | [Dataset tools](dataset_tools.md) |
| `benchmark/` | Time the model and compare kernel variants | [Benchmarks](benchmarks.md) |
| `benchmark/profile.jl` | Find where the time and memory of a run go | [Profiling](profiling.md) |
| `benchmark/fortran/` | Check the model against the original Fortran and time both | [Fortran comparison](fortran.md) |
| `tools/diagnostics/` | Measure a property of the model, such as its budgets | [Diagnostics](diagnostics.md) |
| `tools/validation/` | Check that a change leaves the results unchanged | [Validation](validation.md) |
| `test/` | Run the test suite and write new tests | [Testing](testing.md) |
| `test/support/preset_reference.jl` | Regenerate the preset reference after a change made on purpose | [Testing](testing.md) |
| `.claude/skills/docs-check/check_docs.jl` | Check the local links and version claims in the Markdown docs | [Code conventions](conventions.md) |
| - | Float32, allocation, `@turbo`, formatting and the docs build | [Code conventions](conventions.md) |

Experiments that answer a climate question, such as the effective climate sensitivity, are in the separate repository [GREBExperiments](https://github.com/EnvDroneSense/GREBExperiments), not here.

Scripts that run the model read the local dataset with `greb_data_dir(; allow_download=false)` and never download it.
