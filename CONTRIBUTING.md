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
│   │   ├── calendar.jl         # ModelTime, step_of_year, decimal_year
│   │   ├── postprocess.jl      # monthly climatology/anomalies
│   │   ├── model.jl            # init_model!/qflux_correction!/greb_model!
│   │   └── ensemble.jl         # run_ensemble
│   ├── config/                 # Config parts, presets, RunSpec, resolve()
│   ├── forcing/                # forcing(): CO2 and sunlight per timestep
│   ├── physics/                # radiation.jl, hydrology.jl, ocean.jl, circulation.jl
│   ├── data.jl                 # greb_data_dir(): dataset location + DataDep
│   └── io.jl                   # JLD2 loaders
├── ext/                        # GREBEnzymeExt.jl: Enzyme rules, loaded only with Enzyme
├── test/                       # one file per subject; SHARD in support/selection.jl lists them
│   └── data/                   # golden and preset references (see below)
├── benchmark/                  # timing harness (run_benchmarks.jl, README.md)
├── docs/                       # Documenter site
├── examples/                   # user examples: quickstart.jl, feedback_deconstruction.jl
├── viz/                        # plotting toolbox + Pluto notebook (notebooks/), own Project.toml
├── tools/                      # maintainer scripts: dataset/, diagnostics/, experiments/, validation/ (see tools/README.md)
├── DATA_README.md              # raw .bin input inventory (maintainers)
└── CHANGELOG.md
```

`greb_input_data/` (the dataset) and `Data/` (raw `.bin` inputs) are expected
at runtime or by maintainers but gitignored.

## Input data

The package needs the GREB input dataset (439 MB unpacked, a 321 MB download) to run the model.
`greb_data_dir()` resolves it in this order:

1. an explicit path you pass
2. the `GREB_DATA` environment variable
3. `greb_input_data/` in the package's top folder
4. a copy downloaded earlier
5. a download via [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl)

Only step 5 touches the network, and it prompts before downloading. See
[DATA_README.md](DATA_README.md) for the raw-file inventory.

You do **not** need the dataset to contribute. Tests that require it resolve
with `allow_download=false` and skip when it is absent, which is exactly how
CI runs.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

While developing, run only what a change needs: `julia --project=. test/support/changed.jl --run`. The whole suite is the gate before you open a pull request.

- A new test file needs a row in the `SHARD` table in [test/support/selection.jl](test/support/selection.jl), or it will not run.
- A test that calls `greb_data_dir` must pass `allow_download=false`.
- CI has no dataset, so the golden, MSCM, deconstruction and orbital-table tests skip there. Run the whole suite locally with the dataset.

The developer docs page [Testing](docs/src/dev/testing.md) has every command and the rules for new tests; [test/README.md](test/README.md) maps the files.

## Conventions

Fields are `Float32`, kernels do not allocate, and four kernels stay free of `@turbo`. A change that must not alter results is checked with `tools/validation/bit_identity.jl`. Source folders follow `.JuliaFormatter.toml`, and a CI job fails on unformatted files. See [Code conventions](docs/src/dev/conventions.md) for the details, the formatting command and the docs build.

## Branches

| Branch | Role |
|---|---|
| `main` | Released code only |
| `dev` | Integration; fixes that are not physics modules collect here |
| `dev-refactor` | Structural (non-physics) refactoring, from `dev` |
| `refactor/<item>` | Optional item branch, from `dev-refactor`; merges back into it |
| `physics/<module>` | One physics module each, from `dev` |
| `archive/*` (tags) | Old branch tips kept for reference |

Open pull requests against the branch you branched from. Shared
infrastructure (hooks, loaders, formats) goes on `dev-refactor` and reaches
`dev` before a module relies on it; keep infrastructure and module code in
separate commits. A `physics/*` branch can be frozen when nothing else
depends on it; `physics/volcanic` is frozen, so do not build on it.

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
  [Model overview](https://EnvDroneSense.github.io/GREBClimate.jl/stable/model/).

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
