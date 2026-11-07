---
name: benchmark
description: Run GREBClimate.jl's timing benchmarks (benchmark/run_benchmarks.jl) and its profiler (benchmark/profile.jl) and report results honestly, accounting for machine noise. Use when the user asks to benchmark, time, profile, or check the speed of the model, or to compare thread counts. Allocations are checked by the test suite, not here.
---

# GREB benchmark

`benchmark/run_benchmarks.jl` is a harness with five modes, all on the real dataset after JIT warm-up. `benchmark/README.md` lists the modes, flags, defaults and set-up line; this file holds the procedure. Benchmarks are not a correctness signal; the test suite is, and allocations are checked by `test/test_invariants.jl`.

| Mode | Answers | Command |
|---|---|---|
| `year` (default) | How fast is a 1-year `:full_model` control run; before vs after a change | `julia --project=benchmark -t 2,0 benchmark/run_benchmarks.jl year` |
| `stages` | Where the time goes in six timed calls (not a whole step) | `julia --project=benchmark -t 2,0 benchmark/run_benchmarks.jl stages` |
| `threads` | Speedup of `-t 1..4`, each in its own subprocess | `julia --project=benchmark -t 2,0 benchmark/run_benchmarks.jl threads` |
| `years` | Cost and stability over many years: `--ctrl=N --scnr=N [--experiment=NAME]` | `julia --project=benchmark -t 2,0 benchmark/run_benchmarks.jl years --ctrl=10 --scnr=100` |
| `compare` | Same-process, interleaved comparison with a self-check (`compare_kernels` in `benchmark/compare.jl` takes any variants) | `taskset -c <core> julia --project=benchmark -t 2,0 benchmark/run_benchmarks.jl compare` |

An unknown mode, a mistyped flag or a missing dataset stops the run with an error.

## Reading a run

The harness prints the machine state itself and marks a run NOISY; a NOISY run is not quoted. Each run appends a CSV line per variant to `benchmark/results/` (gitignored). Numbers from different sessions are not comparable (see `.claude/rules/benchmarks.md`). Check the power line first: on battery a `year` reads several times slower than on mains, and a faster reading is not a speedup if the load was simply lighter.

## Steps

1. Pick the mode that answers the question (table above). Use `-t 2,0` (no interactive thread; the `threads` mode launches `-t N,0` too).
2. Before trusting a slow or surprising `year`/`threads` reading, rule out noise:

   | Cause | Check |
   |---|---|
   | Stale precompile cache | `julia --project -e 'using Pkg; Pkg.precompile()'`; if it recompiles `GREBClimate`, re-run |
   | Post-reboot background load | `tasklist \| grep -i searchindexer`; wait 1-2 minutes |
   | Normal wall-clock swing | 50-150 ms across repeated `year` runs; average several runs or sweeps |

   `stages` is far less noise-prone; a single reading there is more trustworthy.
3. For thread counts use `threads`, and run it several times: the threaded paths carry most of the variance (`-t 1` is stable). Compare relative speedups, not absolute times.
4. Report min/median (`year`, `threads`), or the per-stage table (`stages`). Say plainly when noise makes a reading untrustworthy.

## Profiling

`benchmark/profile.jl` shows where the time goes; the harness above shows how much there is. Reports land in `benchmark/profiles/<date>-<commit>-<mode>/` (gitignored).

| Mode | Answers | Command |
|---|---|---|
| `step` (default) | Which function and line of the time loop owns the time | `julia --project=. -t 1 benchmark/profile.jl step` |
| `setup` | What one `greb_model!` call costs besides the time loop, in ms | `julia --project=. -t 1 benchmark/profile.jl setup` |
| `allocs` | Which lines allocate, and how many bytes | `julia --project=. -t 1 benchmark/profile.jl allocs` |
| `dispatch` | Whether a call in the step is dispatched at run time (static; needs JET.jl in the default environment) | `julia --project=. benchmark/profile.jl dispatch` |
| `compare` | What changed in the shares between two `step` or `setup` runs | `julia --project=. benchmark/profile.jl compare <before> <after>` |

| Rule | Why |
|---|---|
| Run at `-t 1` | On Windows the sampler records the first thread only |
| Check `header.txt` first: 10000 samples or more | The sampler's rate varies (2 to 13 ms per sample); the script repeats the run up to a target and warns below 5000 |
| Read `category.txt`, then `owned.txt`, then `tree.txt` | Kind of cost, then the owning line, then the call path |
| Treat a row below three standard errors (`+/-`) as noise | Also for the `diff` column of `compare` |
| A `@turbo` loop is one row at its `@turbo for` line | Lines inside it are not resolved |
| In `compare`, read the row that rose | The other rows fall by dilution without having changed |
| A share is a ceiling, not a gain | Prove a gain with the same-process comparison at the end of this file |

Reference shape of `step` (2026-10-05, release 2.0.1, `-t 1`): `circulation!` 89.5% on the stack (`_diffusion!` 40%, `_advection!` 35% self), 87% of samples inside `@turbo` code, no runtime dispatch, allocation only in `output!`.

## Other checks

| Check | When |
|---|---|
| `InteractiveUtils.code_native`, grep for `gather` and `cvtss2sd`/`cvtsd2ss` | A change touches a `@turbo` loop; a width mismatch or stray `Float64` can hide behind a wall-clock number |
| `test/test_threading.jl` "threaded circulation matches serial" | After touching the parallel branch of `tendencies!` or either `circulation!` call; `Pkg.test()` alone is single-threaded |
| `test/test_golden.jl` or a diff against a saved run | With every timing change; a faster wrong answer is not a win |

To compare code variants use `compare_kernels` (it shuffles the order, times a second copy of the baseline and prints `NOT OK` when that copy is more than 1 % off 1.000; discard such a run). `--counters` adds instructions and cycles per call on Linux and needs LinuxPerf in the user environment, `perf_event_paranoid` <= 2 and a pinned P-core; it refuses a value whose counter did not run the whole time. On Windows there are no counters and no pinning, so effects of a few percent at year scale cannot be resolved.
