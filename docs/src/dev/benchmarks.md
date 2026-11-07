# Benchmarks

`benchmark/run_benchmarks.jl` times the model. It asserts nothing about results: the test suite is the correctness gate, and allocations are checked by `test/test_invariants.jl`. The full list of modes and flags is in `benchmark/README.md`.

| Mode | Measures |
|:-----|:---------|
| `year` (default) | One control year of `:full_model` on the stored flux corrections, 3 runs |
| `stages` | One call of each of six physics stages, about 0.5 s each |
| `threads` | `year` in a fresh process at `-t 1` to `4` |
| `years` | `--ctrl=N` control and `--scnr=N` scenario years of `--experiment=<preset>` |
| `compare` | Two variants of a kernel in one process, interleaved, with a self-check |

## In the REPL

Start Julia in the `benchmark` environment (set up once, see `benchmark/README.md`), load the file and call `main`. Loading defines `main` and runs no mode. From a shell the same run is `julia --project=benchmark -t 2,0 benchmark/run_benchmarks.jl year`.

```julia
julia> include("benchmark/run_benchmarks.jl")

julia> main(["year"])
Threads.nthreads() = 2
ctrl=1 scnr=0 (1 simulated years/rep), experiment=full_model
machine: processor load 56 %, mains, background: OneDrive,SearchIndexer
calibration loop: 40.47 ms, 1.61x the best seen here (25.07 ms)
  run 1: 0.498 s  (0.498 s/simulated year)
  run 2: 0.482 s  (0.482 s/simulated year)
  run 3: 0.533 s  (0.533 s/simulated year)
min: 0.482 s  median: 0.498 s  = 0.498 s/simulated year (median) over 1 years
median / minimum of the runs: 1.03
calibration loop after: 29.64 ms, 1.18x the best
machine state: NOISY (calibration 1.61x the best) - do not quote these timings
0
```

`main` returns the exit code, which is the `0` at the end. The figures above come from a busy machine and are shown for the layout, not as results.

## Reading the output

| Line | Meaning |
|:-----|:--------|
| `calibration loop` | A fixed loop, timed before and after the run and compared with the best time this machine has given for it. A ratio above 1.05 marks the run NOISY |
| `median / minimum of the runs` | Above 1.05 also marks the run NOISY |
| `machine state: quiet` or `NOISY (...)` | Quote the timings only when quiet; the reason is in brackets |
| `self-check ... OK` | `compare` only: the copy of the baseline measured within 1 % of 1.000 |

Every run also appends one line per variant to `benchmark/results/<date>.csv` (git-ignored).

## Per stage and comparison

`stages` prints the minimum and median of each stage and its share of the sum of the minima:

```julia
julia> main(["stages"])
Per-stage timing (about 0.5 s each, single workspace, 2 thread(s) available but unused here):
  circulation!(Ta)  min  180.50 µs   median  203.25 µs   ( 48.3% of measured total)
  circulation!(q)   min  173.00 µs   median  189.35 µs   ( 46.3% of measured total)
  SWradiation!      min    2.28 µs   median    2.40 µs   (  0.6% of measured total)
  LWradiation!      min   10.25 µs   median   10.90 µs   (  2.7% of measured total)
  hydro!            min    5.55 µs   median    5.90 µs   (  1.5% of measured total)
  deep_ocean!       min    1.88 µs   median    2.20 µs   (  0.5% of measured total)
  measured total (sum of stage minima): 373.4 µs
```

`compare` times the variants in a shuffled order each round and a second copy of the baseline as a control. It exits with 1 when that control is more than 1 % off:

```julia
julia> main(["compare", "--seconds=2"])
variant                      min (us)   x base   med (us)   x base
circulation!(Ta)              160.300    1.000    168.100    1.000
circulation!(Ta) (copy)       160.400    1.001    168.050    1.000
circulation!(q)               168.600    1.052    176.500    1.050
self-check (copy / baseline, median; min gives 1.001): 1.0  OK  [20 samples each]
```

For other variants call `compare_kernels` (`benchmark/compare.jl`) from a script.

## Errors

An unknown mode, a flag the mode does not take, a reps argument where none is allowed and a missing dataset stop the run, with exit code 1 from a shell:

```julia
julia> main(["year", "--secnds=1"])
ERROR: the year mode takes no --secnds flag
```

## Rules

| Rule | Why |
|:-----|:----|
| Use `-t 2,0` | Two compute threads and no interactive thread |
| Do not compare timings between sessions | An untouched kernel has measured 0.45 us in one session and 1.28 us in another |
| Show a speed-up only with `compare_kernels` | One process, shuffled order, and a control copy of the baseline |
| Do not record test-suite timings as benchmark results | Wall-clock times are not stable |
