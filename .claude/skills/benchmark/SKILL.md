---
name: benchmark
description: Run GREBClimate.jl's timing/allocation benchmarks (benchmark/run_benchmarks.jl) and report results honestly, accounting for machine noise. Use when the user asks to benchmark, time, profile, or check performance/allocations of the model, or to compare thread counts.
---

# GREB benchmark

`benchmark/run_benchmarks.jl` is a dependency-free harness with five modes, all on the real dataset after JIT warm-up. Benchmarks are not a correctness signal; the test suite is.

| Mode | Answers | Command |
|---|---|---|
| `year` (default) | How fast is a 1-year `:full_model` control run; before vs after a change | `julia --project=. -t 2 benchmark/run_benchmarks.jl year` |
| `stages` | Where the time goes: each `tendencies!` stage and its share | `julia --project=. -t 1 benchmark/run_benchmarks.jl stages` |
| `threads` | Speedup of `-t 1..4`, each in its own subprocess | `julia --project=. -t 1 benchmark/run_benchmarks.jl threads` |
| `alloc` | Bytes allocated by one `tendencies!` call; budget 256 bytes (`TENDENCIES_ALLOC_BUDGET`, same value in `test/test_invariants.jl`) | `julia --project=. -t 1 benchmark/run_benchmarks.jl alloc` |
| `years` | Cost and stability over many years: `--ctrl=N --scnr=N [--experiment=NAME]` (default 10/10, `full_model`, which holds CO2 flat; use `a1b_scenario` for a trend) | `julia --project=. -t 2 benchmark/run_benchmarks.jl years --ctrl=10 --scnr=100` |

The data directory is the next positional argument (default `../greb_input_data`, or `GREB_DATA`); `[jld2_dir] [reps]` are positional in every mode. The legacy form `run_benchmarks.jl <dir>` runs `year`.

## Reference numbers (this machine)

| Regime | `year` at `-t 2` | Recorded |
|---|---|---|
| Normal background load | ~0.6-0.75 s per simulated year | 2026-08-21, mean 0.63 s |
| Background processes closed | ~0.27-0.43 s per simulated year | 2026-09-22, mean 0.31 s |

| Measure | Value (2026-09-22) |
|---|---|
| `circulation!` share of a timestep | 94.5% (Ta and q about half each) |
| `-t 2` over `-t 1` | 1.13x mean (1.01-1.26x); still the best count |
| `-t 3`, `-t 4` | 0.87-1.00x; never help on this grid |
| Allocations per `tendencies!` | 0 bytes |

Check which regime the machine is in (`tasklist`) before judging a `year` reading: 0.6 s is normal under load, and 0.3 s is not a speedup if the load was simply lighter.

## Steps

1. Pick the mode that answers the question (table above). Default to `-t 2`.
2. Before trusting a slow or surprising `year`/`threads` reading, rule out noise:

   | Cause | Check |
   |---|---|
   | Stale precompile cache | `julia --project -e 'using Pkg; Pkg.precompile()'`; if it recompiles `GREBClimate`, re-run |
   | Post-reboot background load | `tasklist \| grep -i searchindexer`; wait 1-2 minutes |
   | Normal wall-clock swing | 50-150 ms across repeated `year` runs; average several runs or sweeps |

   `stages` and `alloc` are far less noise-prone; a single reading there is more trustworthy.
3. For thread counts use `threads`, and run it several times: the threaded paths carry most of the variance (`-t 1` is stable). Compare relative speedups, not absolute times.
4. Report mean/min/max (`year`, `threads`), the per-stage table (`stages`) or the byte count (`alloc`). Say plainly when noise makes a reading untrustworthy.

## Other checks

| Check | When |
|---|---|
| `InteractiveUtils.code_native`, grep for `gather` and `cvtss2sd`/`cvtsd2ss` | A change touches a `@turbo` loop; a width mismatch or stray `Float64` can hide behind a wall-clock number |
| Isolated `@benchmark` of one function (BenchmarkTools in a scratch environment) with realistic input | A change scoped to one kernel; confirm the magnitude afterwards with `year`/`stages` |
| `test/test_threading.jl` "threaded circulation matches serial" | After touching the parallel branch of `tendencies!` or either `circulation!` call; `Pkg.test()` alone is single-threaded |
| `test/test_golden.jl` or a diff against a saved run | With every timing change; a faster wrong answer is not a win |

For measurements beyond this harness: compile old and new into one process, alternate them in shuffled order, take the minimum of N trials, and time a second copy of the baseline as a self-check (it must land at about 1.00x).
