---
paths: ["benchmark/**/*.jl"]
---

`benchmark/README.md` documents every mode, the noise rules and the profiler.
The rules that matter while editing here:

- **Wall-clock times from different sessions are not comparable** on this
  machine (an untouched `convergence!` read 0.45 µs in one session and 1.28 µs
  in another). A speed-up is shown only with `compare_kernels` in
  `benchmark/compare.jl`: one process, shuffled order each round, and a second
  copy of the baseline as a self-check. Discard a run whose self-check fails.
- Linux can pin (`taskset -c <P-core>`) and count (`--counters`); Windows can do
  neither, so effects below about 5 % of a year cannot be shown there.
- Time with `-t 2,0`. Profile with `-t 1`: on Windows the sampler records the
  first thread only.
- A profile share is the ceiling of a gain, not a gain. A row smaller than
  three of its standard errors (the `+/-` column) is noise.
- Chairmarks lives only in `benchmark/Project.toml`. `common.jl` must stay free
  of it, because `profile.jl` includes it under `--project=.`.
- Never record test-suite wall-clock times as benchmark results.
