---
paths: ["benchmark/**/*.jl"]
---

**Cross-session wall-clock comparison on this machine is worthless** — an
untouched `convergence!` measured 0.45 µs in one session and 1.28 µs in
another with zero code change. Compile old and new kernels into one process
and interleave the trials, shuffling order within each trial (a fixed order
measurably penalizes whichever variant runs last). Always time a second
identical copy of the baseline as a control; if that self-check is not
~1.00x, the run is not trustworthy and gets thrown away, not reported.

Use `-t 2,0` for the `year`/`stages` modes (two compute threads, no
interactive thread: 0.259 s per year against 0.303 s for plain `-t 2`,
measured 2026-10-01). Re-measured 2026-09-22 on a quiet machine, plain `-t 2`
was the best count but only 1.01-1.26x over `-t 1` (mean 1.13x);
`-t 3`/`-t 4` never help. The gain may depend on background load, so judge a
threading change with a same-process comparison.

Never record wall-clock test timings as if they were benchmark results.
Assertion counts from the test suite are stable across runs; timings from
this harness are not — don't conflate the two kinds of number.

To find where the time goes, use `benchmark/profile.jl` (modes `step`,
`setup`, `allocs`, `dispatch`, `compare`; see `benchmark/README.md`) and read
`category.txt`, then `owned.txt`. Run it at `-t 1`: on Windows the sampler
records the first thread only. A share in a profile is the ceiling of a gain,
not a measurement of one - a row smaller than three of its standard errors
(the `+/-` column) is noise, and a speed-up still needs the same-process
comparison above.
