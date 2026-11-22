---
paths: ["benchmark/**/*.jl"]
---

`benchmark/README.md` documents every mode, the noise rules and the profiler.
The rules that matter while editing here:

- **Wall-clock times from different sessions are not comparable** on this
  machine (an untouched `convergence!` read 0.45 µs in one session and 1.28 µs
  in another). A speed-up is shown only with `compare` (`benchmark/run.jl compare`
  or `Bench.compare`): one process, shuffled order each round, and a second
  copy of the baseline as a self-check. Discard a run whose self-check fails.
- Linux can pin (`taskset -c <P-core>`) and count (`--counters`); Windows can do
  neither, so effects below about 5 % of a year cannot be shown there.
- Time with `-t 2,0`. Profile with `-t 1`: on Windows the sampler records the
  first thread only.
- A profile share is the ceiling of a gain, not a gain. A row smaller than
  three of its standard errors (the `+/-` column) is noise.
- Chairmarks lives only in `benchmark/Project.toml`. `tools/ToolBox.jl` must stay
  free of it, and the `driver` mode of the Fortran script loads the package only,
  so the Julia time it reports has no benchmark packages in it.
- Never record test-suite wall-clock times as benchmark results.
