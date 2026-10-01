# Contributing to GREBClimate.jl

Contributions are welcome - bug reports, physics questions, performance work,
documentation fixes. This page covers how to get the package running locally,
what the test suite expects, and the few conventions that are easy to trip over.

## Getting set up

```bash
git clone https://github.com/EnvDroneSense/GREBClimate.jl.git
cd GREBClimate.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Project structure

```
GREBClimate.jl/
├── src/                        # the package (module GREBClimate)
│   ├── GREBClimate.jl          # module shell + include order
│   ├── core/                   # the model's backbone
│   │   ├── constants.jl        # grid/physical constants
│   │   ├── state.jl            # ClimateFields, ModelState, workspaces
│   │   ├── tendencies.jl       # per-timestep physics pipeline
│   │   ├── budgets.jl          # shared flux sums, BudgetCheck (the observer's budget check)
│   │   ├── output.jl           # diagnostics!/output!/time_loop!
│   │   ├── postprocess.jl      # monthly climatology/anomalies
│   │   └── model.jl            # init_model!/qflux_correction!/greb_model!
│   ├── config/                 # Config parts, presets, RunSpec, resolve()
│   ├── forcing/                # forcing(): CO2 and sunlight per timestep
│   ├── physics/                # radiation.jl, hydrology.jl, ocean.jl, circulation.jl
│   ├── data.jl                 # greb_data_dir(): dataset location + DataDep
│   └── io.jl                   # JLD2 loaders
├── test/                       # one file per subject; runtests.jl lists them
│   └── data/                   # generated preset reference (see below)
├── benchmark/                  # timing/allocation suite (run_benchmarks.jl, README.md)
├── docs/                       # Documenter site
├── examples/                   # plain-Julia drivers: run_greb.jl, parameter_sweep.jl
├── notebooks/                  # Pluto explorer + launcher (uses the viz/ environment)
├── viz/                        # plotting toolbox, own Project.toml
├── tools/                      # maintainer scripts: dataset/, diagnostics/, validation/ (see tools/README.md)
├── DATA_README.md              # raw .bin input inventory (maintainers)
└── CHANGELOG.md
```

`greb_input_data/` (the dataset) and `Data/` (raw `.bin` inputs) are expected
at runtime or by maintainers but gitignored.

## Input data

The package needs the GREB input climatology (about 353 MB) to run the model.
`greb_data_dir()` resolves it in this order:

1. an explicit path you pass
2. the `GREB_DATA` environment variable
3. a local `greb_input_data/` directory
4. a download via [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl)

Only step 4 touches the network, and it prompts before downloading. See
[DATA_README.md](DATA_README.md) for the raw-file inventory.

You do **not** need the dataset to contribute. Tests that require it resolve
with `allow_download=false` and skip when it is absent, which is exactly how
CI runs.

## Running the tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The suite is split into two shards, grouped by measured runtime rather than
file count and only roughly balanced. CI runs them as separate jobs; locally
you can run one:

```bash
GREB_TEST_SHARD=light julia --project=. -e 'using Pkg; Pkg.test()'
GREB_TEST_SHARD=heavy julia --project=. -e 'using Pkg; Pkg.test()'
```

Things worth knowing before you add tests:

- **A new test file must be added to the `SHARD` table in
  [test/runtests.jl](test/runtests.jl), or it will not run.** Nothing globs
  the directory.
- Shared fixtures live in `test/testutils.jl`: `quiet()`, `with_tempdir()`,
  `synthetic_fields()`, `constant_fields()`, `uniform_record()`,
  `at_first_step()`, `DATA_DIR`, and the grid constants.
- Any test that calls `greb_data_dir` must pass `allow_download=false`.
  Without it the test passes locally (where a dataset usually exists) and
  fails in CI.
- **Don't run the model to test something the model doesn't do.** Most of the
  suite's cost is simulated years. If the behaviour under test lives in
  `init_model!`, a loader, or a single physics kernel, call that directly.
  If it is decided at the start of a run (the CO2 or solar table a preset
  uses, the climatology in use), `at_first_step()` runs `greb_model!` to the
  first step through the observer, returns what it saw and stops. Pass
  `NoCorrections()` there: the observer is not called during the spin-up.
- A test should fail when the behaviour breaks. Before deleting or adding a
  test, break the line it claims to pin and check that it fails (and that, for
  a deletion, another test does too).
- `test/test_golden.jl` is a regression lock on model output. If your change
  moves those numbers, that is either a bug or a deliberate physics change
  that needs saying out loud in the PR.
- **CI has no dataset, so the golden, MSCM, deconstruction and orbital-table
  tests skip there.** A green CI says nothing about them: with the dataset
  present, run the full suite locally before opening a PR.
- `test/test_presets.jl` compares every preset (its forcing, masks, physics,
  hydrology and corrections) with `test/data/preset_reference.jl`. If you
  change a preset on purpose, regenerate the file in the same commit with
  `julia --project=. tools/validation/preset_reference.jl` and say so in the
  PR; never edit it by hand.

## Performance conventions

The model is written to run many simulated years, and the kernels are built
around that:

- Fields are `Float32` throughout, on a fixed `xdim x ydim` grid.
- Physics kernels write into pre-allocated `CirculationWorkspace` buffers
  instead of allocating. `test/test_invariants.jl` enforces a small byte
  budget per physics kernel and checks its return type is concrete - a change
  that allocates per grid cell will fail it immediately. A new kernel must be
  added to both tables in that file. Only the allocation table is checked
  for completeness (an assertion fails if an exported `!` kernel is missing);
  the return-type table is not, so add the kernel there yourself.
- Configuration is passed explicitly. Nothing is held as module-global
  mutable state.

**A change that must not alter results** (a refactor, an optimisation) is
checked with exact equality, not tolerances. With the dataset present, save a
snapshot on the commit before and compare after, at one and at two threads:

```bash
julia --project=. -t 1 tools/validation/bit_identity.jl save snap_t1.jld2
julia --project=. -t 1 tools/validation/bit_identity.jl compare snap_t1.jld2
```

Snapshots depend on the CPU and the Julia version and are not committed. For
a change to the flux sums or the humidity limiters,
`tools/diagnostics/budget.jl <preset>` prints how well each store closes.

**Benchmarks:**

```bash
julia --project=. -t 2 benchmark/run_benchmarks.jl year
```

Two threads is the current recommendation: the temperature and humidity
transport run in parallel, and a third thread has measured no consistent gain.
It has not been re-measured since the circulation was last optimised.

## Documentation

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Docstring examples in ```` ```jldoctest ```` blocks are executed during the
docs build, and a mismatch fails the build. If you change an exported
function's behaviour or signature, check whether an example needs updating.

## Branches

| Branch | Role |
|---|---|
| `main` | Released code only |
| `dev` | Integration; fixes that are not physics modules collect here |
| `dev-refactor` | Structural (non-physics) refactoring, from `dev` |
| `physics/<module>` | One physics module each, from `dev` |

Open pull requests against the branch you branched from. Shared
infrastructure (hooks, loaders, formats) goes on `dev-refactor` and reaches
`dev` before a module relies on it; keep infrastructure and module code in
separate commits.

## Opening a pull request

- Say what changed and why; if the physics moved, say which numbers moved.
- Keep the test suite green, including the golden regression.
- New behaviour needs a test. New exported names need a docstring.
- Note user-visible changes in [CHANGELOG.md](CHANGELOG.md).

## Reporting a bug

First check that the input data is complete and where `greb_data_dir()`
expects it, and that your Julia and package versions meet the requirements.
Then open an issue with the Julia version (`versioninfo()`), the OS, the
experiment configuration (`preset`/`Config` call and any options you
changed), and the full error or the unexpected output. A `RunSpec` short enough
to reproduce quickly helps a lot.

## Roadmap

Contributions towards these are especially welcome:

- **NetCDF output** - optional direct write of monthly means.
- **Visualisation dashboard** - interactive maps and time series, similar to the
  [MSCM interactive database](https://mscm.dkrz.de/GREB_model.html?locale=EN).
- **Physics guide** - a derivation-level companion to the
  [Model overview](https://EnvDroneSense.github.io/GREBClimate.jl/dev/model/).

## Credits

GREBClimate.jl is a Julia translation. The model itself is the work of the
GREB developers at Monash University:

- **Dietmar Dommenget** - original GREB model
- **Janine Flöter** - original GREB model
- **Tobias Bayr** - GREB development
- **Christian Stassen** - hydrological cycle (MSCM)
- **Kerry Nice, Mike Rezny, Dietmar Kasang** - Monash Simple Climate Model
  experiments and database

See the References section of the [README](README.md) for the papers, which are
the right thing to cite for the model.

The Julia package is the work of:

- **Thomas Struys** (UGent) - Julia translation and optimization
- **Michiel Stock** (UGent) - Julia development guidance, initial package refactor

If you contribute, add yourself here.
