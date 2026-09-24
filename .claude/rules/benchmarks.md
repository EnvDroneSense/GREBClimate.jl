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

Use `-t 2` for the `year`/`stages` modes. Re-measured 2026-09-22 on a quiet
machine: `-t 2` is still the best count but only 1.01-1.26x over `-t 1`
(mean 1.13x); `-t 3`/`-t 4` never help. The gain may depend on background
load, so judge a threading change with a same-process comparison.

Never record wall-clock test timings as if they were benchmark results.
Assertion counts from the test suite are stable across runs; timings from
this harness are not — don't conflate the two kinds of number.
