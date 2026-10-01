---
paths: ["test/**/*.jl"]
---

`test/test_invariants.jl` builds its coverage guard (the `all_kernels` set) from
`names(GREBClimate)` — only *exported* `!`-kernels are ever picked up. An
unexported kernel (e.g. `_diffusion!`, `_advection!`) is invisible to the
guard directly; it is only covered because the exported wrapper that calls it
(`diffusion!`, `advection!`) is in the `kernels` allocation table. A new
unexported kernel with no exported caller would silently never be checked —
this asymmetry is easy to get wrong.

The `kernels` table is guarded by an equality assertion (`setdiff(all_kernels,
...)`); the return-type `signatures` table has **no** such assertion — adding a
new kernel there is convention, not machine-checked.

`Pkg.test()` forces `--check-bounds=yes` (see the comment on `exe` in
`test_threading.jl`).
This is deliberate, not overhead to shave off: the code indexes ghost-cell
buffers by hand under `@inbounds`/`@turbo`, exactly where bounds checking
earns its cost. Do not disable or "optimize" it away.

Shared fixtures live in `test/testutils.jl` (`synthetic_fields`,
`constant_fields`, `uniform_record`, `quiet`, `with_tempdir`, `gmean`).
`at_first_step` runs `greb_model!` to the first step of a phase through the
observer and returns what it saw there: use it for what a preset decides at
the start (the CO2, the solar table, the climatology in use) instead of
simulating a year; it is not called during the spin-up, so pass
`NoCorrections()` or `SpinUp(0)`. Prefer exact `==` over `isapprox` where the claim
is bit-identity — see the `-t 1` vs `-t 2` monthly-mean comparison in
`test_threading.jl`. Do not run the model to test something the model does
not do.
