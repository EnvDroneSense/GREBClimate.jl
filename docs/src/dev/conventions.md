# Code conventions

The model runs many simulated years, so its kernels are written around that. These are the rules that are easy to break by accident.

| Convention | Detail |
|:-----------|:-------|
| `Float32` throughout | Fields are `Float32` on the fixed `xdim x ydim` grid. A stray `Float64` literal slows a kernel down |
| No allocation in a kernel | Physics kernels write into pre-allocated `ModelWorkspace` buffers. `test/test_invariants.jl` fails a kernel that allocates per grid cell |
| No `@turbo` in four kernels | `LWradiation!`, `hydro!`, `SWradiation!` and `seaice!` use plain loops so that Enzyme can differentiate them. Their `log` and `exp` go through `fastlog` and `fastexp`, whose derivative rules are in `ext/GREBEnzymeExt.jl`. `circulation!` and the other kernels use `@turbo` |
| Explicit configuration | Configuration is passed as an argument. Nothing is held as module-global mutable state |

## A change that must not alter results

A refactor or an optimisation is checked with exact equality, not with tolerances. With the dataset present, save a snapshot on the commit before the change and compare after it, once at one thread and once at two:

```bash
julia --project=. -t 1 tools/validation/bit_identity.jl save snap_t1.jld2
julia --project=. -t 1 tools/validation/bit_identity.jl compare snap_t1.jld2
```

Repeat both with `-t 2,0` and a second snapshot. Snapshots depend on the CPU and the Julia version, so they are not committed. For a change to the flux sums or the humidity limiters, `tools/diagnostics/budget.jl <preset>` prints how well each store closes. See [Validation](validation.md) and [Diagnostics](diagnostics.md).

## Speed

[Benchmarks](benchmarks.md) time the model, [Profiling](profiling.md) shows where the time goes, and [Fortran comparison](fortran.md) checks the model against the original. Two threads is the current recommendation: the temperature and humidity transport run in parallel, and a third thread has measured no consistent gain. It has not been re-measured since the circulation was last optimised.

## Formatting

`src/`, `test/`, `benchmark/`, `tools/`, `ext/`, `examples/` and `viz/` follow `.JuliaFormatter.toml`, and a CI job fails when a file is not formatted. With JuliaFormatter 2 in your global environment, format before you push:

```bash
julia -e 'using JuliaFormatter; foreach(format, ("src", "test", "benchmark", "tools", "ext", "examples", "viz"))'
```

## Documentation

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Docstring examples in ```` ```jldoctest ```` blocks run during the docs build, and a mismatch fails the build. If you change the behaviour or signature of an exported function, check whether an example needs updating. `julia .claude/skills/docs-check/check_docs.jl` checks links, version claims and exported names in seconds.

The site has two live versions: [stable](https://EnvDroneSense.github.io/GREBClimate.jl/stable/) is the latest release, built when the release tag is created, and [dev](https://EnvDroneSense.github.io/GREBClimate.jl/dev/) follows the `dev` branch. Pull requests build the docs and run the doctests but deploy nothing.
