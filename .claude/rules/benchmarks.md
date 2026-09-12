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

Use `-t 2` for the `year`/`stages` modes. This was tuned when `circulation!`
was ~98% of per-timestep cost; it is now ~93% after the ghost-cell change
and has not been re-measured since — treat `-t 2` as the current
recommendation, not a settled conclusion.

Never record wall-clock test timings as if they were benchmark results.
Assertion counts from the test suite are stable across runs; timings from
this harness are not — don't conflate the two kinds of number.
