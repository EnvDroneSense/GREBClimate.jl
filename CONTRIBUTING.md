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
  `synthetic_fields()`, `DATA_DIR`, and the grid constants.
- Any test that calls `greb_data_dir` must pass `allow_download=false`.
  Without it the test passes locally (where a dataset usually exists) and
  fails in CI.
- **Don't run the model to test something the model doesn't do.** Most of the
  suite's cost is simulated years. If the behaviour under test lives in
  `init_model!`, a loader, or a single physics kernel, call that directly.
- `test/test_golden.jl` is a regression lock on model output. If your change
  moves those numbers, that is either a bug or a deliberate physics change
  that needs saying out loud in the PR.

## Performance conventions

The model is written to run many simulated years, and the kernels are built
around that:

- Fields are `Float32` throughout, on a fixed `xdim x ydim` grid.
- Physics kernels write into pre-allocated `CirculationWorkspace` buffers
  instead of allocating. `test/test_invariants.jl` enforces a small byte
  budget per physics kernel and checks its return type is concrete - a change
  that allocates per grid cell will fail it immediately. A new kernel must be
  added to both tables in that file; an assertion there fails if you forget.
- Configuration is passed explicitly. Nothing is held as module-global
  mutable state.

Benchmarks:

```bash
julia --project=. -t 2 benchmark/run_benchmarks.jl year
```

Two threads, not three: circulation was ~98% of per-timestep cost pre-ghost-cell;
after `865ae01`'s periodic ghost cells it is ~93% and has not been re-measured
since. Either way there is no third lane of reliable work.

## Documentation

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Docstring examples in ```` ```jldoctest ```` blocks are executed during the
docs build, and a mismatch fails the build. If you change an exported
function's behaviour or signature, check whether an example needs updating.

## Opening a pull request

- Say what changed and why; if the physics moved, say which numbers moved.
- Keep the test suite green, including the golden regression.
- New behaviour needs a test. New exported names need a docstring.
- Note user-visible changes in [CHANGELOG.md](CHANGELOG.md).

## Reporting a bug

Open an issue with the Julia version, the OS, the experiment configuration
(`create_experiment_config` call and any switches you changed), and the full
error or the unexpected output. A `RunSpec` short enough to reproduce quickly
helps a lot.

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
