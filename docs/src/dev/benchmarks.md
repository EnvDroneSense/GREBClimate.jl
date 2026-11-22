# Benchmarks

`benchmark/run.jl` times the model. It asserts nothing about results: the test suite is the correctness gate, and allocations are checked by `test/test_invariants.jl`. The full list of flags is in `benchmark/README.md`. To find which function or line holds the time, see [Profiling](profiling.md). To compare with the original Fortran model, see [Fortran comparison](fortran.md).

```bash
julia --project=benchmark -t 2,0 benchmark/run.jl <mode> [--flag=value ...]
```

| Mode | Measures |
|:-----|:---------|
| `time` (default) | `--ctrl=N` control and `--scnr=N` scenario years, `--reps=N` times (1, 0 and 3) |
| `threads` | The same run in a fresh process at `-t 1` to `4` |
| `kernels` | Each kernel of one step, about 0.5 s each, and `step!` as a whole |
| `compare` | Two or more kernels (`--only=A,B`) in one process, interleaved, with a self-check |

Every mode takes `--preset=<name>` (default `full_model`). The profiling modes are on the next page.

## In the REPL

Start Julia in the `benchmark` environment (set up once, see `benchmark/README.md`), load the file and call `main`. Loading defines `main` and runs no mode.

```julia
julia> include("benchmark/run.jl")

julia> main(["time"])
machine: processor load 21 %, mains, background: OneDrive,SearchIndexer
calibration loop: 26.76 ms, 1.07x the best seen here (25.07 ms)
  run 1: 0.507 s  (0.507 s/simulated year)
  run 2: 0.549 s  (0.549 s/simulated year)
  run 3: 0.512 s  (0.512 s/simulated year)
min: 0.507 s  median: 0.512 s  = 0.512 s/simulated year (median) over 1 years
median / minimum of the runs: 1.01
calibration loop after: 30.87 ms, 1.23x the best
machine state: NOISY (calibration 1.23x the best) - do not quote these timings
0
```

`main` returns the exit code, which is the `0` at the end. The figures come from a busy machine and are shown for the layout, not as results.

## Reading the output

| Line | Meaning |
|:-----|:--------|
| `calibration loop` | A fixed loop, timed before and after the run and compared with the best time this machine has given for it. A ratio above 1.05 marks the run NOISY |
| `median / minimum of the runs` | Above 1.05 also marks the run NOISY |
| `machine state: quiet` or `NOISY (...)` | Quote the timings only when quiet; the reason is in brackets |
| `self-check ... OK` | `compare` only: the copy of the baseline measured within 1 % of 1.000 |

Every run also appends one line per variant to `benchmark/results/<date>.csv` (git-ignored).

## Per kernel and comparison

`kernels` prints the minimum and median of each kernel and its share of `step!`, one whole step:

```julia
julia> main(["kernels"])
Kernel timing (about 0.5 s each, single workspace, 2 thread(s) available but unused here):
  circulation!(Ta)  min  159.20 µs   median  160.10 µs   ( 45.3% of step!)
  circulation!(q)   min  165.90 µs   median  167.05 µs   ( 47.2% of step!)
  SWradiation!      min    1.90 µs   median    1.91 µs   (  0.5% of step!)
  LWradiation!      min    8.60 µs   median    8.80 µs   (  2.4% of step!)
  hydro!            min    4.60 µs   median    4.68 µs   (  1.3% of step!)
  deep_ocean!       min    1.50 µs   median    1.56 µs   (  0.4% of step!)
  seaice!           min    2.30 µs   median    2.37 µs   (  0.7% of step!)
  step!             min  351.50 µs   median  412.20 µs   (100.0% of step!)
  (convergence! is inside circulation!(q); step! also holds the update loop, output! and
   diagnostics!. For shares of a whole run use the profile mode)
```

`compare` times the variants in a shuffled order each round and a second copy of the baseline as a control. It exits with 1 when that control is more than 1 % off:

```julia
julia> main(["compare", "--only=circulation!(Ta),circulation!(q)", "--seconds=2"])
variant                      min (us)   x base   med (us)   x base
circulation!(Ta)              163.000    1.000    163.450    1.000
circulation!(Ta) (copy)       163.100    1.001    163.750    1.002
circulation!(q)               171.500    1.052    171.950    1.052
self-check (copy / baseline, median; min gives 1.001): 1.002  OK  [20 samples each]
```

For other variants call `Bench.compare` with closures (`benchmark/Bench.jl`):

```julia
Bench.compare("old" => () -> f_old(x), "new" => () -> f_new(x); seconds=5)
```

## Errors

An unknown mode, a flag the mode does not take and a missing dataset stop the run, with exit code 1 from a shell:

```julia
julia> main(["time", "--secnds=1"])
ERROR: the time mode takes no --secnds flag
```

## Rules

| Rule | Why |
|:-----|:----|
| Use `-t 2,0` | Two compute threads and no interactive thread |
| Do not compare timings between sessions | An untouched kernel has measured 0.45 us in one session and 1.28 us in another |
| Show a speed-up only with `compare` | One process, shuffled order, and a control copy of the baseline |
| Do not record test-suite timings as benchmark results | Wall-clock times are not stable |
