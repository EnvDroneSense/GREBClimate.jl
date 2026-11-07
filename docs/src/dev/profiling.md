# Profiling

`benchmark/profile.jl` shows where the time and the memory of a run go. [Benchmarks](benchmarks.md) say how much time there is; the profiler says which function and line holds it. It runs in the package environment and needs no extra packages, except for the two optional modes below. The full list of flags is in `benchmark/README.md`.

| Mode | Measures |
|:-----|:---------|
| `step` (default) | Sampled control runs of `:full_model` (25 simulated years each, repeated until 10000 samples): where the time loop spends its time |
| `setup` | The same with 1-year runs: what one `greb_model!` call costs besides the time loop |
| `allocs` | Every allocation of one 1-year run, by function and line |
| `dispatch` | Whether a call in one model step is dispatched at run time (static; needs JET.jl in the default environment) |
| `compare` | The change in shares between two `step` or `setup` runs |

## In the REPL

Start Julia with one thread in the package environment, load the file and call `main`. Loading defines `main` and runs no mode. This file and `run_benchmarks.jl` both define `main`, so use one per session. From a shell the same run is `julia --project=. -t 1 benchmark/profile.jl step --years=1 --samples=500`.

```julia
julia> include("benchmark/profile.jl")

julia> main(["step", "--years=1", "--samples=500", "--out=run-a"])
# Profile: step

| Item    | Value                                         |
|:--------|:----------------------------------------------|
| mode    | step (full_model control, stored corrections) |
| commit  | e493d38 (tracked files modified)              |
| date    | 2026-11-07 18:39                              |
| julia   | 1.12.6                                        |
| threads | 1 (sampled: 1)                                |
| runs    | 2 of 1 simulated year(s)                      |
| elapsed | 1.51 s (757.2 ms per run)                     |
| samples | 675 (2.24 ms per sample)                      |

## By class

Each sample is classed by its innermost frame. Shares are of 675 samples; +/- is one standard error in percentage points.

| Class                           | Samples | Share % |  +/- |
|:--------------------------------|--------:|--------:|-----:|
| vectorized kernel               |     564 |   83.56 | 1.43 |
| Base (array access, scalar ops) |     101 |   14.96 | 1.37 |
| package code                    |       4 |    0.59 | 0.30 |
| runtime and system (C)          |       3 |    0.44 | 0.26 |
| copy and fill                   |       2 |    0.30 | 0.21 |
| scalar math                     |       1 |    0.15 | 0.15 |

## By function

Self: the function holds the innermost package frame. On stack: it is anywhere in the stack. ms: the on-stack share of one run (757.2 ms).

| Function          | Self | Share % |  +/- | On stack | Share % |     ms |
|:------------------|-----:|--------:|-----:|---------:|--------:|-------:|
| `_diffusion!`     |  293 |   43.41 | 1.91 |      298 |   44.15 | 334.28 |
| `_advection!`     |  196 |   29.04 | 1.75 |      205 |   30.37 | 229.96 |
| `circulation!`    |   67 |    9.93 | 1.15 |      592 |   87.70 | 664.07 |
| `_convergence!`   |   17 |    2.52 | 0.60 |       17 |    2.52 |  19.07 |
| `refresh_ghosts!` |   16 |    2.37 | 0.59 |       16 |    2.37 |  17.95 |
| `hydro!`          |   11 |    1.63 | 0.49 |       14 |    2.07 |  15.70 |
| `LWradiation!`    |   10 |    1.48 | 0.47 |       20 |    2.96 |  22.43 |
| `fastlog`         |   10 |    1.48 | 0.47 |       10 |    1.48 |  11.22 |

## By line

The innermost package frame. A `@turbo` loop is one row, at its `@turbo for` line.

| Line                               | Self | Share % |  +/- |
|:-----------------------------------|-----:|--------:|-----:|
| `circulation.jl:178  _diffusion!`  |  113 |   16.74 | 1.44 |
| `circulation.jl:140  _diffusion!`  |   99 |   14.67 | 1.36 |
| `circulation.jl:279  _advection!`  |   86 |   12.74 | 1.28 |
| `circulation.jl:433  circulation!` |   59 |    8.74 | 1.09 |
| `circulation.jl:129  _diffusion!`  |   50 |    7.41 | 1.01 |
| `circulation.jl:315  _advection!`  |   45 |    6.67 | 0.96 |
| `circulation.jl:344  _advection!`  |   33 |    4.89 | 0.83 |
| `circulation.jl:61  _convergence!` |   16 |    2.37 | 0.59 |

full report: run-a\report.md
0
```

`main` returns the exit code, which is the `0` at the end. This run is short (675 samples, where the default takes 10000) and is shown for the layout, not as a result: the profiler warns below 5000 samples. The screen shows the top eight rows of each table; `report.md` in the run's folder has all of them.

## Reading the report

| Section | Read it for |
|:--------|:------------|
| By class | What kind of time there is: vectorized kernel, Base array access, package code, allocation and GC, copies, dispatch, compilation, waiting, file I/O |
| By function | Which function owns the time. A sample belongs to the innermost frame in `src/`, so library time is charged to the line that called it |
| By line | Which line of the model owns the time |
| Largest frames (`report.md` only) | Under the class table: the frames of each class and the function each is charged to |

The other files of a `step` or `setup` run:

| File | Content |
|:-----|:--------|
| `counts.tsv` | Self samples of every class, function and line, uncut, as `table`, `row`, `value` columns; `compare` reads it |
| `tree.txt` | The call tree from `greb_model!` down, rows below 0.1 percent left out |
| `flat.txt` | Self time per frame, libraries included |
| `stacks.folded` | One line per distinct stack, for flame graph tools (speedscope, FlameGraph, inferno) |
| `profile.pb.gz` | The same samples in pprof format, with `--viewer=pprof` (needs PProf.jl in the default environment) |

Reports go to `benchmark/profiles/<date>-<commit>-<mode>/`, which is git-ignored, or to `--out`.

## Comparing two runs

```julia
julia> main(["compare", "run-a", "run-b"])
# Profile comparison

| Run | Directory | Commit                           | Samples |
|:----|:----------|:---------------------------------|--------:|
| A   | run-a     | e493d38 (tracked files modified) |     675 |
| B   | run-b     | e493d38 (tracked files modified) |     501 |
...
## By function

| Row             |   A % |   B % |  Diff |  +/- | 3 SE |
|:----------------|------:|------:|------:|-----:|-----:|
| `circulation!`  |  9.93 |  7.19 | -2.74 | 1.63 |      |
| `_advection!`   | 29.04 | 31.14 | +2.10 | 2.71 |      |
...
```

A run is a report directory or its name under `benchmark/profiles/`. Rows are sorted by the size of the change; a `*` in the last column marks a difference of three standard errors or more. A row one run lacks counts as zero there. A run made before `counts.tsv` existed must be profiled again.

## Errors

An unknown mode, a flag the mode does not take, a second positional argument, a missing dataset and a missing `counts.tsv` stop the run, with exit code 1 from a shell:

```julia
julia> main(["step", "--yeers=1"])
ERROR: the step mode takes no --yeers flag
```

## Rules

| Rule | Why |
|:-----|:----|
| Run at `-t 1` | On Windows the sampler records the first thread only, so at `-t 2,0` one of the two `circulation!` calls is missing |
| Take 5000 samples or more | The sampler's rate varies from 2 to 13 ms per sample; the script repeats the run up to the target and warns below 5000 |
| A row smaller than three standard errors is noise | The same holds for the Diff column of `compare` |
| A share is the ceiling of a gain, not a gain | A profile shows where the time goes. Show that a change is faster with `compare_kernels` (see [Benchmarks](benchmarks.md)) |
| When one row grows, the other shares fall | In `compare`, read the row that rose and ignore the rest |
| A `@turbo` loop is one row | The whole loop is charged to its `@turbo for` line; lines inside it are not resolved |
