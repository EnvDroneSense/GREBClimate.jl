# Benchmarks

Timing and allocation harness for GREBClimate.jl. It asserts nothing about
results: the test suite is the correctness gate, and a faster benchmark with a
red suite is a regression. It needs the local dataset and never downloads it
(`GREB_DATA` or a path argument overrides the location).

The timing modes need Chairmarks, which lives only in `benchmark/Project.toml`
and is not a dependency of the package. Set the environment up once:

```bash
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=benchmark -t 2,0 benchmark/run_benchmarks.jl [mode] [jld2_dir] [reps]
```

| Mode | Measures | Default reps | Notes |
|:-----|:---------|:-------------|:------|
| `year` (default) | One control year of `:full_model` on the stored flux corrections | 3 | Prints each run and the min and median in seconds. `--out=PATH` also writes the run times to a file, one per line |
| `stages` | One call of each physics stage: `circulation!` (air, vapor), `SWradiation!`, `LWradiation!`, `hydro!`, `deep_ocean!` | none | `--seconds=S` per stage (default 0.5) instead of a rep count. Prints min and median per call; the share is of the minima. Single workspace, so threads are not used. Leaves out `seaice!`, output and the tendency assembly |
| `threads` | `year` in a fresh process at `-t 1`, 2, 3, 4 | 3 | Prints the speedup against `-t 1` |
| `alloc` | Bytes allocated by one `tendencies!` call | none | Compared with the 256-byte budget in `test/test_invariants.jl`; the two numbers are kept in sync by hand |
| `years` | `--ctrl=N` control years plus `--scnr=N` scenario years of `--experiment=<preset>` | 1 | Defaults 10, 10, `full_model`. Stored corrections, no spin-up. Prints seconds per simulated year |
| `compare` | `circulation!(Ta)` against `circulation!(q)` in one process, interleaved | none | `--seconds=S` per variant (default 5), `--counters` (below). A fixed demo of `compare_kernels`; exits 1 when the self-check fails |

Example: `julia --project=benchmark -t 2,0 benchmark/run_benchmarks.jl years --ctrl=10 --scnr=100 --experiment=co2_double`.

## Reading the numbers

| Rule | Why |
|:-----|:----|
| Use `-t 2,0` for `year`, `years` and `stages` | Two compute threads and no interactive thread. Measured 2026-10-01 at 0.259 s per year against 0.303 s for plain `-t 2` and 0.407 s for `-t 1`; `-t 3` and `-t 4` do not help. CI uses the same |
| Do not compare timings between sessions | An untouched kernel has measured 0.45 us in one session and 1.28 us in another |
| Compare variants in one process | Compile both, interleave the trials, shuffle the order each trial, and time a second copy of the baseline as a control. If the control is not about 1.00x, discard the run |
| Do not record test-suite timings as benchmark results | Assertion counts are stable; wall-clock times are not |

`year`, `years`, `stages` and `compare` print the machine state before and
after: the processor load, the power source (on Linux the power profile, the
CPU governor and the CPU the process runs on), and a fixed calibration loop compared with
the best time this machine has given for it (kept in
`benchmark/profiles/calibration.txt`, which is gitignored). A run is marked
NOISY when the calibration loop is more than 5 % off the best, or when the
median of the runs is more than 5 % above their minimum, or when the Linux
power profile is `power-saver`; do not quote a noisy run.

`stages` gives shares of its six calls, not of a step. For the shares of a
whole run, including what `stages` leaves out, use the `step` profile below.

`year`, `years` and `stages` time one variant per run. They are a quick check
and a source of per-stage shares, not a measurement protocol: a claimed
speed-up needs the same-process comparison above, which is `compare`.

## Comparing variants

`compare_kernels` in `compare.jl` implements that comparison:

```julia
include("benchmark/common.jl"); include("benchmark/compare.jl")
compare_kernels("old" => () -> f_old(x), "new" => () -> f_new(x); seconds=5)
```

The first pair is the baseline. A second copy of it is timed as a control. The
run is `rounds` (20) rounds; in each, every variant is timed for
`seconds / rounds` in a freshly shuffled order. A sample is the fastest call of
one round's chunk, so the table shows the min and median of 20 samples and the
ratio of each to the baseline. If the copy's median ratio is not within 1 % of
1.000 the function prints `NOT OK` and returns `false`; discard that run. Each
call compiles its own variants first, so define them in the same session.

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
| Pinning | Pin to a P-core: `taskset -c <P-core> julia --project=benchmark -t 1 benchmark/run_benchmarks.jl compare --counters` |

When a requirement is missing, the flag prints one line saying which and the
run goes on without counters.

On Windows there are no counters and no pinning, so timings are noisier;
effects below about 5 % at the scale of a year cannot be resolved there.

## Results

Every run of `year`, `years`, `stages` and `compare` appends a line per variant
to `benchmark/results/<date>.csv` (gitignored): `date,commit,mode,name,threads,samples,min_s,median_s,ratio_min,ratio_median,selfcheck,state`.
The ratio columns are empty outside `compare`, `selfcheck` is filled on the
baseline-copy row, and `state` is `quiet` or `NOISY`. The commit carries
`-dirty` when tracked files are modified.

## Profiling

```bash
julia --project=. -t 1 benchmark/profile.jl [mode] [jld2_dir] [--years=N] [--samples=N] [--out=DIR] [--viewer=pprof]
```

| Mode | Profiles | Default | Answers |
|:-----|:---------|:--------|:--------|
| `step` (default) | Sampled control runs of `:full_model` on the stored flux corrections, repeated until `--samples` samples are collected (at most 8 runs) | 25 years, 10000 samples | Where the time loop spends its time |
| `setup` | The same with short runs (at most 400) | 1 year, 10000 samples | What one `greb_model!` call costs besides the time loop |
| `allocs` | Every allocation of one run | 1 year | Which lines allocate, and how much |
| `dispatch` | Nothing is run: JET's optimization analysis of one `time_loop!` step | - | Whether any call in the step is dispatched at run time |

Reports go to `benchmark/profiles/<date>-<commit>-<mode>/`, which is gitignored:

| File | Modes | Content | Read it for |
|:-----|:------|:--------|:------------|
| `header.txt` | all | Commit, Julia version, threads, runs, elapsed time, sample count and interval | Whether the run is usable |
| `category.txt` | `step`, `setup` | Samples per kind of cost, judged by the innermost frame: vectorized kernel, Base array access, package code, allocation and GC, copies, dispatch, compilation, waiting, I/O | What kind of time there is |
| `owned.txt` | `step`, `setup` | Samples per package function and per package line; a sample belongs to its innermost frame in `src/`, so library time is charged to the line that called it. The `ms` column is the function's share of one run | Which line of the model owns the time |
| `flat.txt` | `step`, `setup` | Self time per frame, libraries included, largest last | The single hottest frames |
| `tree.txt` | `step`, `setup` | The call tree from `greb_model!` down, rows below 0.1 percent left out | The call path behind a row |
| `stacks.folded` | `step`, `setup` | One line per distinct stack: the frames from `greb_model!` inwards separated by `;`, then the sample count | A flame graph: load the file in speedscope or any tool that reads collapsed stacks. Also plain text to search |
| `profile.pb.gz` | `step`, `setup`, with `--viewer=pprof` | The same samples in pprof format | The pprof web view, and its `-list=<function>` report of counts beside the source |
| `dispatch.txt` | `dispatch` | The number of runtime-dispatch reports inside the package, and each report with its call chain | Type instabilities in the step; zero reports is the expected state |
| `allocs.txt` | `allocs` | Bytes and counts per package function and per package line, with the type holding the most bytes | Where memory is allocated |

| Rule | Why |
|:-----|:----|
| Run it at `-t 1` | On Windows the sampler records the first thread only, so at `-t 2,0` one of the two `circulation!` calls is missing |
| Read the `+/-` column | It is one standard error. A row smaller than three of them, or a difference between two runs smaller than that, is noise |
| A `@turbo` loop is one row | The whole loop is charged to its `@turbo for` line; lines inside it are not resolved |
| The `ms` column is a mean | It includes the garbage collection a function causes, so an allocating function reads higher here than in a best-of-N timing |
| A share is not a speed | A profile shows where the time goes. That a change is faster needs the same-process comparison above |

`--viewer=pprof` needs PProf.jl and the `dispatch` mode needs JET.jl in the
default Julia environment (`julia -e 'using Pkg; Pkg.add("PProf")'`); neither
is a dependency of the package. Open the file with `using PProf; PProf.refresh(file="<path>/profile.pb.gz")`.

To see what changed between two `step` or `setup` runs:

```bash
julia --project=. benchmark/profile.jl compare <run before> <run after>
```

A run is a report directory or its name under `benchmark/profiles/`. The
output lists, by class, function and line, each row's share in both runs, the
difference in percentage points and its standard error, largest change first;
`*` marks a difference of three standard errors or more. When one row grows,
the shares of all others fall, so read the row that rose and ignore the rest.

The sampler's rate varies on this machine (2 to 13 ms per sample has been
seen), which is why a run is repeated up to a sample target. In `allocs.txt`,
`Profile.Allocs.BufferType` is the data buffer of an array.

## Files

| File | Role |
|:-----|:-----|
| `run_benchmarks.jl` | The modes above |
| `compare.jl` | `compare_kernels`: the interleaved comparison |
| `counters.jl` | Hardware counters through LinuxPerf; loaded only by `--counters` |
| `profile.jl` | The sampling profile above |
| `common.jl` | Argument parsing, dataset lookup, machine state and result records, shared with `profile.jl` and the Fortran comparison. Needs no Chairmarks |
| `Project.toml` | Chairmarks and the package for `run_benchmarks.jl`. `Manifest.toml` is gitignored. `profile.jl` still runs in the package environment (`--project=.`) |
| `fortran/` | Local only (excluded in `.git/info/exclude`, not in the repository). `run_fortran_comparison.jl` times and checks GREBClimate.jl against the original Fortran GREB (modes `compare`, `memory`, `verify`, `io`, `build`, `selftest`); `peak_memory.ps1` is its Windows peak-memory helper. Needs gfortran and the Fortran source (`GREB_GFORTRAN`, `GREB_FORTRAN_DIR`) |
