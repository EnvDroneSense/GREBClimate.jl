# Benchmarks

Timing and allocation harness for GREBClimate.jl. It asserts nothing about
results: the test suite is the correctness gate, and a faster benchmark with a
red suite is a regression. It needs the local dataset and never downloads it
(`GREB_DATA` or a path argument overrides the location).

```bash
julia --project=. -t 2,0 benchmark/run_benchmarks.jl [mode] [jld2_dir] [reps]
```

| Mode | Measures | Default reps | Notes |
|:-----|:---------|:-------------|:------|
| `year` (default) | One control year of `:full_model` on the stored flux corrections | 3 | Prints each run and mean/min/max in seconds |
| `stages` | One call of each physics stage: `circulation!` (air, vapour), `SWradiation!`, `LWradiation!`, `hydro!`, `deep_ocean!` | 2000 | Single workspace, so threads are not used. Leaves out `seaice!`, output and the tendency assembly |
| `threads` | `year` in a fresh process at `-t 1`, 2, 3, 4 | 3 | Prints the speedup against `-t 1` |
| `alloc` | Bytes allocated by one `tendencies!` call | none | Compared with the 256-byte budget in `test/test_invariants.jl`; the two numbers are kept in sync by hand |
| `years` | `--ctrl=N` control years plus `--scnr=N` scenario years of `--experiment=<preset>` | 1 | Defaults 10, 10, `full_model`. Stored corrections, no spin-up. Prints seconds per simulated year |

Example: `julia --project=. -t 2,0 benchmark/run_benchmarks.jl years --ctrl=10 --scnr=100 --experiment=co2_double`.

## Reading the numbers

| Rule | Why |
|:-----|:----|
| Use `-t 2,0` for `year`, `years` and `stages` | Two compute threads and no interactive thread. Measured 2026-10-01 at 0.259 s per year against 0.303 s for plain `-t 2` and 0.407 s for `-t 1`; `-t 3` and `-t 4` do not help. CI uses the same |
| Do not compare timings between sessions | An untouched kernel has measured 0.45 us in one session and 1.28 us in another |
| Compare variants in one process | Compile both, interleave the trials, shuffle the order each trial, and time a second copy of the baseline as a control. If the control is not about 1.00x, discard the run |
| Do not record test-suite timings as benchmark results | Assertion counts are stable; wall-clock times are not |

This harness times one variant per run. It is a quick check and a source of
per-stage shares, not a measurement protocol: a claimed speed-up needs the
same-process comparison above.

## Profiling

```bash
julia --project=. -t 1 benchmark/profile.jl step [jld2_dir] [--years=N] [--samples=N] [--out=DIR]
```

`step` samples control runs of `:full_model` on the stored flux corrections
(`--years`, default 50), repeated until `--samples` samples are collected
(default 10000, at most 8 runs), and writes five files to
`benchmark/profiles/<date>-<commit>-step/`, which is gitignored:

| File | Content | Read it for |
|:-----|:--------|:------------|
| `header.txt` | Commit, Julia version, threads, years, elapsed time, sample count and interval | Whether the run is usable |
| `category.txt` | Samples per kind of cost, judged by the innermost frame: vectorized kernel, Base array access, package code, allocation and GC, copies, dispatch, compilation, waiting, I/O | What kind of time there is |
| `owned.txt` | Samples per package function and per package line; a sample belongs to its innermost frame in `src/`, so library time is charged to the line that called it | Which line of the model owns the time |
| `flat.txt` | Self time per frame, libraries included, largest last | The single hottest frames |
| `tree.txt` | The call tree from `greb_model!` down, rows below 0.1 percent left out | The call path behind a row |

| Rule | Why |
|:-----|:----|
| Run it at `-t 1` | On Windows the sampler records the first thread only, so at `-t 2,0` one of the two `circulation!` calls is missing |
| Read the `+/-` column | It is one standard error. A row smaller than three of them, or a difference between two runs smaller than that, is noise |
| A `@turbo` loop is one row | The whole loop is charged to its `@turbo for` line; lines inside it are not resolved |
| A share is not a speed | A profile shows where the time goes. That a change is faster needs the same-process comparison above |

The sampler's rate varies on this machine (2 to 13 ms per sample has been
seen), which is why the run is repeated up to a sample target.

## Files

| File | Role |
|:-----|:-----|
| `run_benchmarks.jl` | The modes above |
| `profile.jl` | The sampling profile above |
| `common.jl` | Argument parsing and dataset lookup shared with the Fortran comparison |
| `Manifest.toml` | Gitignored; there is no `Project.toml` here, scripts run in the package environment (`--project=.` from the repo root) |
| `fortran/` | Local only (excluded in `.git/info/exclude`, not in the repository). `run_fortran_comparison.jl` times and checks GREBClimate.jl against the original Fortran GREB (modes `compare`, `memory`, `verify`, `io`, `build`, `selftest`); `peak_memory.ps1` is its Windows peak-memory helper. Needs gfortran and the Fortran source (`GREB_GFORTRAN`, `GREB_FORTRAN_DIR`) |
