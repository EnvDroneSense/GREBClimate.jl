# Benchmarks

Timing, comparison and profiling harness for GREBClimate.jl. It asserts nothing about
results: the test suite is the correctness gate, and a faster benchmark with a
red suite is a regression. It needs the local dataset and never downloads it
(`GREB_DATA` or `--data=DIR` overrides the location).

It needs Chairmarks, which lives only in `benchmark/Project.toml` and is not a
dependency of the package. Set the environment up once, and again after updating
the repository (it is not rebuilt for you):

```bash
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=benchmark -t 2,0 benchmark/run.jl <mode> [--flag=value ...]
```

| Mode | Does | Flags (default) |
|:-----|:-----|:----------------|
| `time` (default) | Times control (and scenario) runs; prints each run, the min and median, and seconds per simulated year | `--ctrl` (1), `--scnr` (0), `--reps` (3), `--out=FILE` (run times, one per line) |
| `threads` | The same run at `-t 1`, 2, 3, 4, each in its own process; prints the speedup against `-t 1` | `--threads=1,2`, `--ctrl`, `--scnr`, `--reps` |
| `ab` | The `time` mode of an older commit against this checkout, in pairs of fresh processes in alternating order; the older commit is a temporary git worktree. Prints the ratio of the medians and a self-check; exits 1 when the base processes differ by more than 5 % | `--base=REV` (required), `--rounds` (3), `--ctrl`, `--scnr`, `--reps` per process |
| `kernels` | Times each kernel of one step, and `step!` (a whole step), with the share of `step!` | `--only=A,B` (all), `--seconds` per kernel (0.5) |
| `compare` | Times two or more kernels interleaved, with a self-check; the first is the baseline. Exits 1 when the self-check fails | `--only=A,B[,C]` (required), `--seconds` per variant (5), `--counters` |
| `profile` | Samples control runs | `--years` (25), `--samples` (10000), `--out=DIR` |
| `allocs` | Records every allocation of one run | `--years` (1), `--out=DIR` |
| `dispatch` | Checks one model step for runtime dispatch with JET (static) | `--out=DIR` |
| `diff` | The change in shares between two profile directories | `diff <before> <after>` |

Every mode but `diff` takes `--preset=NAME` (default `full_model`) and `--data=DIR`. A flag a mode does not take, an unknown mode and a missing dataset stop the run with an error. The kernels are the entries of `Bench.kernel_table`: `circulation!(Ta)`, `circulation!(q)`, `SWradiation!`, `LWradiation!`, `hydro!`, `deep_ocean!`, `seaice!` and `step!`.

Examples:

```bash
julia --project=benchmark -t 2,0 benchmark/run.jl time --ctrl=10 --scnr=100 --preset=co2_double
julia --project=benchmark -t 2,0 benchmark/run.jl compare --only="circulation!(Ta),circulation!(q)"
julia --project=benchmark -t 1 benchmark/run.jl profile
```

Allocations are not timed here: `test/test_invariants.jl` checks every kernel, `tendencies!` included, and runs in CI (`Pkg.test(test_args=["invariants"])`).

## Your own code

`Bench.compare`, `Bench.profile`, `Bench.allocs` and `Bench.dispatch` take closures, so any code can be timed or profiled:

```julia
include("benchmark/Bench.jl")
Bench.compare("old" => () -> f_old(x), "new" => () -> f_new(x); seconds=5)
```

`Bench.profile(f; warmup=false)` skips the warm-up call when `f` is already compiled. Load `Bench.jl` or `run.jl` in a session; `main(["kernels"])` in `run.jl` returns the exit code.

## Reading the numbers

| Rule | Why |
|:-----|:----|
| Use `-t 2,0` for `time`, `kernels` and `compare` | Two compute threads and no interactive thread. Measured 2026-10-01 at 0.259 s per year against 0.303 s for plain `-t 2` and 0.407 s for `-t 1`; `-t 3` and `-t 4` do not help. CI uses the same |
| Do not compare timings between sessions | An untouched kernel has measured 0.45 us in one session and 1.28 us in another |
| Compare variants in one process | `compare` compiles both, interleaves the trials, shuffles the order each round, and times a second copy of the baseline as a control. If the control is not about 1.00x, discard the run |
| Do not record test-suite timings as benchmark results | Assertion counts are stable; wall-clock times are not |

`time`, `kernels`, `compare` and the Fortran comparison print the machine state before and
after: the processor load, the power source (on Linux the power profile, the
CPU governor and the CPU the process runs on), and a fixed calibration loop compared with
the best time this machine has given for it (kept in
`benchmark/results/calibration.txt`, which is gitignored). A run is marked
NOISY when the calibration loop is more than 5 % off the best, or when the
median of the runs is more than 5 % above their minimum, or when the Linux
power profile is `power-saver`; do not quote a noisy run.

`kernels` gives shares of one step, with a single workspace and no threads. For the shares of a
whole run, use `profile`. `time` and `kernels` time one variant per run: a claimed
speed-up needs `compare`.

## Comparing variants

`compare` runs `rounds` (20) rounds. In each, every variant is timed for
`seconds / rounds` in a freshly shuffled order. A sample is the fastest call of
one round's chunk, so the table shows the min and median of 20 samples and the
ratio of each to the baseline. A second copy of the baseline is the control: if its
median ratio is not within 1 % of 1.000 the run prints `NOT OK` and returns `false`
(exit code 1 from a shell); discard that run. Each call compiles its own variants first.

Ratios below about 2 % need `--counters` or a quiet machine. The timings
cannot tell a 1 % change from the machine.

## Counters (Linux only)

`compare --counters` also prints instructions and cycles per call from the
hardware counters. They resolve differences far below what the clock does:
`circulation!(q)` against `circulation!(Ta)` reads 1.030 on instructions and
about 1.047 on cycles.

| Requirement | Detail |
|:------------|:-------|
| LinuxPerf | `julia -e 'using Pkg; Pkg.add("LinuxPerf")'` in your own default environment. It is not in `benchmark/Project.toml`, so Windows resolves cleanly |
| `perf_event_paranoid` | `/proc/sys/kernel/perf_event_paranoid` must be 2 or lower |
| Fill rate | A variant whose counters did not run the whole time (fill rate below 1.000, for example on an E-core) is printed as refused |
| Pinning | Pin to a P-core: `taskset -c <P-core> julia --project=benchmark -t 1 benchmark/run.jl compare --only=... --counters` |

When a requirement is missing, the flag prints one line saying which and the
run goes on without counters.

On Windows there are no counters and no pinning, so timings are noisier;
effects below about 5 % at the scale of a year cannot be resolved there.

## Results

Every run of `time`, `threads`, `ab`, `compare` and the Fortran comparison appends a line per variant
to `benchmark/results/<date>.csv` (gitignored): `date,commit,mode,name,threads,samples,min_s,median_s,ratio_min,ratio_median,selfcheck,state`.
The ratio columns are empty outside `compare`, `selfcheck` is filled on the
baseline-copy row, and `state` is `quiet` or `NOISY`. The commit carries
`-dirty` when tracked files are modified.

## Profiling

Run `profile`, `allocs` and `dispatch` as in the table above, with `-t 1` for `profile`. `profile` samples runs of `--years` control years until `--samples` samples are collected; `allocs` records one run; `dispatch` runs nothing and needs JET.jl in the default Julia environment (`julia -e 'using Pkg; Pkg.add("JET")'`; it is not a dependency of the package).

Reports go to `benchmark/profiles/<date>-<commit>-<mode>/` (or `--out`), which is gitignored. Every mode writes `report.md`: a Markdown file with a table of the run (commit, Julia version, preset, threads, runs, elapsed time, samples) followed by the result tables, padded so it also reads as plain text. The run prints the top eight rows of the same tables.

| Mode | `report.md` holds | Read it for |
|:-----|:------------------|:------------|
| `profile` | Samples by class (the kind of cost: vectorized kernel, Base array access, package code, allocation and GC, copies, dispatch, compilation, waiting, I/O), with the largest frames of each; by package function (self and on-stack samples, and the `ms` of one run); by package line. A sample belongs to its innermost frame in `src/`, so library time is charged to the line that called it | What kind of time there is, and which line of the model owns it |
| `allocs` | Bytes and counts per package function and per package line (the 40 largest), with the type holding the most bytes | Where memory is allocated |
| `dispatch` | The number of runtime-dispatch reports inside the package, and each report with its call chain | Type instabilities in the step; zero reports is the expected state |

`profile` also writes:

| File | Content | Read it for |
|:-----|:--------|:------------|
| `counts.tsv` | Self samples of every class, function and line, uncut, as `table`, `row`, `value` columns, plus the commit and sample count | What `diff` reads; also easy to load into a script |
| `tree.txt` | The call tree from `greb_model!` down, rows below 0.1 percent left out | The call path behind a row |
| `stacks.folded` | One line per distinct stack: the frames from `greb_model!` inwards separated by `;`, then the sample count | A flame graph: load the file in speedscope or any tool that reads collapsed stacks. Also plain text to search |

| Rule | Why |
|:-----|:----|
| Run it at `-t 1` | On Windows the sampler records the first thread only, so at `-t 2,0` one of the two `circulation!` calls is missing |
| Read the `+/-` column | It is one standard error. A row smaller than three of them, or a difference between two runs smaller than that, is noise |
| A `@turbo` loop is one row | The whole loop is charged to its `@turbo for` line; lines inside it are not resolved |
| The `ms` column is a mean | It includes the garbage collection a function causes, so an allocating function reads higher here than in a best-of-N timing |
| A share is not a speed | A profile shows where the time goes. That a change is faster needs `compare` |

To see what changed between two `profile` runs:

```bash
julia --project=benchmark benchmark/run.jl diff <run before> <run after>
```

A run is a report directory. The output has one table each for classes, functions and lines, with each row's
share in both runs, the difference in percentage points and its standard
error, largest change first; a `*` in the last column marks a difference of
three standard errors or more. A row one run lacks counts as zero there. When
one row grows, the shares of all others fall, so read the row that rose and
ignore the rest. The page "Profiling" in the developer docs shows the output.

The sampler's rate varies on this machine (2 to 13 ms per sample has been
seen), which is why a run is repeated up to a sample target. In the `allocs`
report, `Profile.Allocs.BufferType` is the data buffer of an array.

## Files

| File | Role |
|:-----|:-----|
| `run.jl` | The command line: the modes above |
| `Bench.jl` | The module behind it: machine state, result records, `compare`, the kernel table, timing and profiling |
| `counters.jl` | Hardware counters through LinuxPerf; included only by `--counters` |
| `Project.toml` | Chairmarks and the package. `Manifest.toml` is gitignored |
| `fortran/` | `run_fortran_comparison.jl` times and checks GREBClimate.jl against the original Fortran GREB (modes `compare`, `verify`, `io`, `build`; `driver` is the Julia run it starts itself). Run it with `julia --project=benchmark benchmark/fortran/run_fortran_comparison.jl [mode] [--reps=N]`. `verify` is the gate and `compare` runs it first. The modes need gfortran and the Fortran source (`GREB_GFORTRAN`, `GREB_FORTRAN_DIR`), and exit 1 when either is missing. `compare` also prints the peak memory of the Julia run (`Sys.maxrss`) |

The shared helpers (flag parsing, dataset lookup) are in `tools/ToolBox.jl`.
