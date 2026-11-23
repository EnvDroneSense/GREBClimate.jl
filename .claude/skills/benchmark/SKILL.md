---
name: benchmark
description: Run GREBClimate.jl's timing, comparison and profiling harness (benchmark/run.jl) and its Fortran comparison (benchmark/fortran/), and report results honestly, accounting for machine noise. Use when the user asks to benchmark, time, profile, or check the speed of the model, or to compare thread counts. Allocations are checked by the test suite, not here.
---

# GREB benchmark

`benchmark/run.jl` is one command line with eight modes, all on the real dataset after JIT warm-up. `benchmark/README.md` lists the modes, flags, defaults and set-up line; this file holds the procedure. Benchmarks are not a correctness signal; the test suite is, and allocations are checked by `test/test_invariants.jl`. An unknown mode, a mistyped flag or a missing dataset stops the run with an error.

| Question | Command (from the repository root) |
|---|---|
| How fast is a control run; before vs after a change | `julia --project=benchmark -t 2,0 benchmark/run.jl time` |
| Where does the time go in one step | `... -t 2,0 benchmark/run.jl kernels` |
| Is variant B faster than A | `... -t 2,0 benchmark/run.jl compare --only=A,B` (kernel names from `kernels`) or `Bench.compare` with closures |
| Speedup of `-t 1..4` | `... -t 2,0 benchmark/run.jl threads` |
| Where does a whole run spend its time | `julia --project=benchmark -t 1 benchmark/run.jl profile` |
| Which lines allocate | `... -t 1 benchmark/run.jl allocs` |
| Is a call in the step dispatched at run time | `... benchmark/run.jl dispatch` (needs JET.jl in the default environment) |
| What changed between two profiles | `... benchmark/run.jl diff <before> <after>` |

Every mode takes `--preset=NAME` and `--data=DIR`.

## Reading a run

The harness prints the machine state itself and marks a run NOISY; a NOISY run is not quoted. Each timing run appends a CSV line per variant to `benchmark/results/` (gitignored). Numbers from different sessions are not comparable (see `.claude/rules/benchmarks.md`). Check the power line first: on battery a `time` run reads several times slower than on mains, and a faster reading is not a speedup if the load was simply lighter.

## Steps

1. Pick the mode that answers the question (table above). Use `-t 2,0` for timing (no interactive thread; `threads` launches `-t N,0` too) and `-t 1` for `profile`.
2. Before trusting a slow or surprising `time` or `threads` reading, rule out noise:

   | Cause | Check |
   |---|---|
   | Stale precompile cache | `julia --project=benchmark -e 'using Pkg; Pkg.precompile()'`; if it recompiles `GREBClimate`, re-run |
   | Post-reboot background load | `tasklist \| grep -i searchindexer`; wait 1-2 minutes |
   | Normal wall-clock swing | 50-150 ms across repeated `time` runs; average several runs or sweeps |

   `kernels` is far less noise-prone; a single reading there is more trustworthy.
3. For thread counts use `threads`, and run it several times: the threaded paths carry most of the variance (`-t 1` is stable). Compare relative speedups, not absolute times.
4. Report min/median (`time`, `threads`), or the per-kernel table (`kernels`). Say plainly when noise makes a reading untrustworthy.

## Profiling

`profile` shows where the time goes; `time` shows how much there is. `benchmark/README.md` has the flags and the files of a run (`report.md` first); the developer docs page "Profiling" shows the output. Reports land in `benchmark/profiles/<date>-<commit>-<mode>/` (gitignored), or in `--out`.

| Rule | Why |
|---|---|
| Run at `-t 1` | On Windows the sampler records the first thread only |
| Check the samples row of the run table first: 10000 or more | The sampler's rate varies (2 to 13 ms per sample); the profiler repeats the run up to a target and warns below 5000 |
| Read the class table, then the function and line tables, then `tree.txt` | Kind of cost, then the owning function and line, then the call path |
| Treat a row below three standard errors (`+/-`) as noise | Also for the `Diff` column of `diff` |
| A `@turbo` loop is one row at its `@turbo for` line | Lines inside it are not resolved |
| In `diff`, read the row that rose | The other rows fall by dilution without having changed |
| A share is a ceiling, not a gain | Prove a gain with `compare` |

Quote shares only from a run of the default size, with its commit and sample count.

## Fortran comparison

`benchmark/fortran/run_fortran_comparison.jl` checks the model against the original Fortran GREB and times both. It needs gfortran and the Fortran source (`GREB_GFORTRAN`, `GREB_FORTRAN_DIR`) and runs in the `benchmark` environment. The developer docs page "Fortran comparison" has the details.

| Mode | Answers | Command |
|---|---|---|
| `verify` | Whether Julia and the Fortran give the same scenario temperatures (the gate; about 1 minute) | `julia --project=benchmark benchmark/fortran/run_fortran_comparison.jl verify` |
| `compare` | Fortran builds and Julia timed in shuffled trials, after the gate passes; also the Julia peak memory | `... compare [--reps=N]` |
| `io`, `build` | The Fortran's output cost, the compile step | `... io` |

| Rule | Why |
|---|---|
| Run `verify` before quoting any Fortran timing | `compare` runs it first and reports no timings when it fails |
| Discard a `compare` run whose self-check (`control / official`) is off by more than 5 % | The two builds are identical, so the machine was noisy |
| The Julia time is a fresh process | Start-up, package load and compilation are inside it, so it is an end-to-end figure, not a kernel timing |
| Exit code 1 means the mode failed | Missing tooling, a failed gate or a failed self-check |

## Other checks

| Check | When |
|---|---|
| `InteractiveUtils.code_native`, grep for `gather` and `cvtss2sd`/`cvtsd2ss` | A change touches a `@turbo` loop; a width mismatch or stray `Float64` can hide behind a wall-clock number |
| `test/test_threading.jl` "threaded circulation matches serial" | After touching the parallel branch of `tendencies!` or either `circulation!` call; `Pkg.test()` alone is single-threaded |
| `test/test_golden.jl` or a diff against a saved run | With every timing change; a faster wrong answer is not a win |

`compare` shuffles the order, times a second copy of the baseline and prints `NOT OK` when that copy is more than 1 % off 1.000; discard such a run. `--counters` adds instructions and cycles per call on Linux and needs LinuxPerf in the user environment, `perf_event_paranoid` <= 2 and a pinned P-core; it refuses a value whose counter did not run the whole time. On Windows there are no counters and no pinning, so effects of a few percent at year scale cannot be resolved.
