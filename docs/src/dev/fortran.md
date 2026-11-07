# Fortran comparison

`benchmark/fortran/run_fortran_comparison.jl` checks GREBClimate.jl against the original Fortran GREB and times both on the same workload (3 flux-correction, 3 control and 50 scenario years, `co2_double`). It builds three gfortran executables, runs them, and runs Julia in a child process. It asserts nothing in the test suite: the Fortran comparison stays in `benchmark/`. The flags of the other harnesses are in [Benchmarks](benchmarks.md).

| Mode | Does |
|:-----|:-----|
| `compare` (default) | Runs `verify`, then times `official`, `native`, `control` and Julia in shuffled trials (`[reps]`, default 3) |
| `verify` | The gate: do both models give the same scenario temperatures (about 1 minute) |
| `memory` | Peak memory of each build and Julia (Windows only, 1 rep) |
| `io` | The cost of the Fortran's output writes |
| `build` | Compiles the three Fortran builds |
| `selftest` | The script's own logic; needs no compiler |
| `driver` | The Julia run, started by the other modes in a child process |

## Set-up

It runs in the package environment, from the repository root:

```bash
julia --project=. benchmark/fortran/run_fortran_comparison.jl [mode] [reps]
```

| Variable | Meaning | Default |
|:---------|:--------|:--------|
| `GREB_GFORTRAN` | Path to gfortran | `gfortran` on `PATH`, then the Scoop MinGW folder |
| `GREB_FORTRAN_DIR` | The Fortran source tree | A folder `greb-official-official` next to the repository, one level up, or in a sibling `data` folder |
| `GREB_FORTRAN_EXP` | Fortran `EXP` number; the Julia preset follows it | `20` (`co2_double`) |
| `GREB_BENCH_DIR` | Build and run tree | `greb_fortran_bench` in the temp folder |
| `GREB_DATA` | Julia input data | The repository dataset |

Missing gfortran or source tree stops every mode except `selftest` with exit code 1 and a line naming what is missing.

## The gate

`verify` runs the Fortran and Julia once, untimed, and compares the monthly scenario Ts. `compare` runs it first and prints no timings when it fails.

```text
Correctness gate - did both models run the same experiment?
  final-year global-mean Ts anomaly: Fortran 2.7581 K, Julia 2.758 K
  Fortran final-year anomaly (K, at least)     2.76     limit 0.1   PASS
  final-year mean, relative difference         2.8e-5   limit 0.01  PASS
  mean |Fortran - Julia| over all months (K)   0.000254 limit 0.005 PASS
  max  |Fortran - Julia| over all months (K)   0.0258   limit 0.2   PASS
  overall -> PASS
```

| Check | Limit | Catches |
|:------|:------|:--------|
| Fortran final-year anomaly | at least 0.1 K | A run with no signal, which any model matches |
| Relative difference of the final-year global mean | 1 % | A different experiment |
| Mean difference over all months | 0.005 K | Drift |
| Largest difference over all months | 0.2 K | A local error |

## Timing

```text
--- trial 1/1 ---
  control  49.55s
  julia    28.3s
  native   28.3s
  official 50.24s

Results (1 reps):
  official min 50.24  median 50.24  max 50.24
  ...
Fortran : Julia (min of trials)
  official / julia = 1.78x
  native   / julia = 1.0x

Self-check (control / official, identical builds): 0.986x  ok
```

The figures come from a noisy machine (the run printed `NOISY`) and show the layout, not results. `compare` prints the machine state before and after, like the other harnesses, and appends one row per variant to `benchmark/results/<date>.csv` with mode `fortran`.

| Build | Flags |
|:------|:------|
| `official` | `-Ofast -ffast-math -funroll-loops -fopenmp` (the original script) |
| `native` | the same plus `-march=native` |
| `control` | a second copy of `official`, the self-check |

## Rules

| Rule | Why |
|:-----|:----|
| Run `verify` before quoting any timing | Two models that disagree have no meaningful speed ratio |
| Discard a run whose self-check is more than 5 % from 1.00 | `control` and `official` are identical builds, so the machine was noisy |
| The Julia time is a fresh process | Start-up, package load and compilation are inside it: an end-to-end figure, not a kernel timing |
| The Fortran also writes about 88 MB | Use `io` to bound what that costs |
| Exit code 1 means the mode failed | Missing tooling, a failed gate or a failed self-check |
