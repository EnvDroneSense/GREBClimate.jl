---
paths: ["test/**/*.jl"]
---

`test/test_invariants.jl` builds its coverage guard (the `all_kernels` set) from
`names(GREBClimate; all=true)`: every `!`-named function defined in a kernel
file, exported or not, must have a row in the `kernels` allocation table, and
the set must not be empty. `_`-prefixed workers (`_diffusion!`, `_advection!`,
`_convergence!`) are skipped; they are covered through the wrapper that calls
them (`diffusion!`, `advection!`, `convergence!`). A new `_` worker with no
wrapper in the table would silently never be checked - give it one.

The `kernels` table is guarded by an equality assertion (`setdiff(all_kernels,
...)`); the return-type `signatures` table has **no** such assertion — adding a
new kernel there is convention, not machine-checked.

`Pkg.test()` forces `--check-bounds=yes` (see the comment on `exe` in
`test_threading.jl`).
This is deliberate, not overhead to shave off: the code indexes ghost-cell
buffers by hand under `@inbounds`/`@turbo`, exactly where bounds checking
earns its cost. Do not disable or "optimize" it away.

The names the tests call besides the exports (the kernels, the run loop, the
calendar, the constants the kernel tests calculate with) are imported once in
`test/support/imports.jl`; add a name there rather than writing
`GREBClimate.name` in many places. `_`-prefixed internals keep the prefix.

Shared fixtures live in `test/support/testutils.jl` (`synthetic_fields`,
`constant_fields`, `uniform_record`, `field`, `quiet`, `with_tempdir`,
`with_dataset`, `gmean`, `run_synthetic`, `record_times`). Input files are
written with `write_field`, `write_flux_corrections`, `write_boundary_anomalies`,
`write_ipcc_scenarios` and `write_solar_scenarios`. One model step is `initial_state(ini)` plus
`step!(t, CO2, s, fields, r)`, which supplies the record buffers and fresh
`state`, `ws` and `acc` unless given.
`run_synthetic(run, config; fields, kwargs...)` is a muted `greb_model!` on
`synthetic_fields()` with `allow_uninitialized`: use it for short model runs
that do not need the dataset.
`at_first_step` runs `greb_model!` to the first step of a phase through the
observer and returns what it saw there: use it for what a preset decides at
the start (the CO2, the solar table, the climatology in use) instead of
simulating a year; it is not called during the spin-up, so pass
`NoCorrections()` or `SpinUp(0)`. Prefer exact `==` over `isapprox` where the claim
is bit-identity — see the `-t 1` vs `-t 2` monthly-mean comparison in
`test_threading.jl`. Do not run the model to test something the model does
not do.

Selecting tests: `SHARD` in `test/support/selection.jl` holds (file, CI shard,
tier, cost in seconds). `test/parallel.jl` uses the cost to balance the groups it
runs side by side; a new file needs a rough measured value. A new test file needs
a row there with the smallest tier that should include it: `smoke` for files under ~10 s, `standard` for kernel, budget and
preset files, `full` for model runs and anything slow. The source-to-test map
is in the root `CLAUDE.md`. A subset run proves only the files it ran.
