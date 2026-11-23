# Profiling

`benchmark/run.jl` shows where the time and the memory of a run go, with the modes `profile`, `allocs`, `dispatch` and `diff`. [Benchmarks](benchmarks.md) say how much time there is; the profiler says which function and line holds it. It runs in the `benchmark` environment. The full list of flags is in `benchmark/README.md`.

| Mode | Measures |
|:-----|:---------|
| `profile` | Sampled control runs of `--preset` (default `full_model`, 25 simulated years each, repeated until 10000 samples): where the time loop spends its time |
| `allocs` | Every allocation of one 1-year run, by function and line |
| `dispatch` | Whether a call in one model step is dispatched at run time (static; needs JET.jl in the default environment) |
| `diff` | The change in shares between two `profile` runs |

## In the REPL

Start Julia with one thread in the `benchmark` environment, load the file and call `main`. Loading defines `main` and runs no mode. From a shell the same run is `julia --project=benchmark -t 1 benchmark/run.jl profile --years=1 --samples=500`.

```julia
julia> include("benchmark/run.jl")

julia> main(["profile", "--years=1", "--samples=500", "--out=run-a"])
# Profile

| Item    | Value                                                                    |
|:--------|:-------------------------------------------------------------------------|
| mode    | profile (full_model, 1 control + 0 scenario year(s), stored corrections) |
| commit  | 6311ed5                                                                  |
| date    | 2026-11-22 22:51                                                         |
| julia   | 1.12.6                                                                   |
| preset  | full_model                                                               |
| years   | 1                                                                        |
| threads | 1 (sampled: 1)                                                           |
| runs    | 3                                                                        |
| elapsed | 1.02 s (338.4 ms per run)                                                |
| samples | 504 (2.01 ms per sample)                                                 |

## By class

Each sample is classed by its innermost frame. Shares are of 504 samples; +/- is one standard error in percentage points.

| Class                           | Samples | Share % |  +/- |
|:--------------------------------|--------:|--------:|-----:|
| vectorized kernel               |     408 |   80.95 | 1.75 |
| Base (array access, scalar ops) |      84 |   16.67 | 1.66 |
| package code                    |       6 |    1.19 | 0.48 |
| copy and fill                   |       5 |    0.99 | 0.44 |
| runtime and system (C)          |       1 |    0.20 | 0.20 |

## By function

Self: the function holds the innermost package frame. On stack: it is anywhere in the stack. ms: the on-stack share of one run (338.4 ms).

| Function          | Self | Share % |  +/- | On stack | Share % |     ms |
|:------------------|-----:|--------:|-----:|---------:|--------:|-------:|
| `_diffusion!`     |  196 |   38.89 | 2.17 |      205 |   40.67 | 137.63 |
| `_advection!`     |  178 |   35.32 | 2.13 |      187 |   37.10 | 125.54 |
| `circulation!`    |   30 |    5.95 | 1.05 |      439 |   87.10 | 294.72 |
| `refresh_ghosts!` |   19 |    3.77 | 0.85 |       19 |    3.77 |  12.76 |
| `_convergence!`   |   12 |    2.38 | 0.68 |       12 |    2.38 |   8.06 |
| `fastlog`         |   11 |    2.18 | 0.65 |       11 |    2.18 |   7.38 |
| `time_loop!`      |    9 |    1.79 | 0.59 |      492 |   97.62 | 330.30 |
| `LWradiation!`    |    9 |    1.79 | 0.59 |       20 |    3.97 |  13.43 |

## By line

The innermost package frame. A `@turbo` loop is one row, at its `@turbo for` line.

| Line                                 | Self | Share % |  +/- |
|:-------------------------------------|-----:|--------:|-----:|
| `circulation.jl:180  _diffusion!`    |   82 |   16.27 | 1.64 |
| `circulation.jl:281  _advection!`    |   64 |   12.70 | 1.48 |
| `circulation.jl:142  _diffusion!`    |   58 |   11.51 | 1.42 |
| `circulation.jl:317  _advection!`    |   40 |    7.94 | 1.20 |
| `circulation.jl:346  _advection!`    |   40 |    7.94 | 1.20 |
| `circulation.jl:131  _diffusion!`    |   29 |    5.75 | 1.04 |
| `circulation.jl:435  circulation!`   |   26 |    5.16 | 0.99 |
| `circulation.jl:17  refresh_ghosts!` |   16 |    3.17 | 0.78 |

full report: run-a\report.md
0
```

`main` returns the exit code, which is the `0` at the end. This run is short (504 samples, where the default takes 10000) and is shown for the layout, not as a result: the profiler warns below 5000 samples. The screen shows the top eight rows of each table; `report.md` in the run's folder has all of them.

## Reading the report

| Section | Read it for |
|:--------|:------------|
| By class | What kind of time there is: vectorized kernel, Base array access, package code, allocation and GC, copies, dispatch, compilation, waiting, file I/O |
| By function | Which function owns the time. A sample belongs to the innermost frame in `src/`, so library time is charged to the line that called it |
| By line | Which line of the model owns the time |
| Largest frames (`report.md` only) | Under the class table: the frames of each class and the function each is charged to |

The other files of a `profile` run:

| File | Content |
|:-----|:--------|
| `counts.tsv` | Self samples of every class, function and line, uncut, as `table`, `row`, `value` columns; `diff` reads it |
| `tree.txt` | The call tree from `greb_model!` down, rows below 0.1 percent left out |
| `stacks.folded` | One line per distinct stack, for flame graph tools (speedscope, FlameGraph, inferno) |

Reports go to `benchmark/profiles/<date>-<commit>-<mode>/`, which is git-ignored, or to `--out`.

## Comparing two runs

```julia
julia> main(["diff", "run-a", "run-b"])
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

A run is a report directory. Rows are sorted by the size of the change; a `*` in the last column marks a difference of three standard errors or more. A row one run lacks counts as zero there.

## Errors

An unknown mode, a flag the mode does not take, a positional argument where none is allowed, a missing dataset and a missing `counts.tsv` stop the run, with exit code 1 from a shell:

```julia
julia> main(["profile", "--yeers=1"])
ERROR: the profile mode takes no --yeers flag
```

## Rules

| Rule | Why |
|:-----|:----|
| Run at `-t 1` | On Windows the sampler records the first thread only, so at `-t 2,0` one of the two `circulation!` calls is missing |
| Take 5000 samples or more | The sampler's rate varies from 2 to 13 ms per sample; the profiler repeats the run up to the target and warns below 5000 |
| A row smaller than three standard errors is noise | The same holds for the Diff column of `diff` |
| A share is the ceiling of a gain, not a gain | A profile shows where the time goes. Show that a change is faster with `compare` (see [Benchmarks](benchmarks.md)) |
| When one row grows, the other shares fall | In `diff`, read the row that rose and ignore the rest |
| A `@turbo` loop is one row | The whole loop is charged to its `@turbo for` line; lines inside it are not resolved |
