---
paths: ["test/**/*.jl"]
---

`test/test_invariants.jl:63-73` builds its coverage guard from
`names(GREBClimate)` — only *exported* `!`-kernels are ever picked up. An
unexported kernel (e.g. `_diffusion!`, `_advection!`) is invisible to the
guard directly; it is only covered because the exported wrapper that calls it
(`diffusion!`, `advection!`) is in the allocation table at `L48-59`. A new
unexported kernel with no exported caller would silently never be checked —
this asymmetry is easy to get wrong.

The allocation table (`L48-59`) is guarded by an equality assertion
(`L73`); the return-type `signatures` table (`L92-102`) has **no** such
assertion — adding a new kernel there is convention, not machine-checked.

`Pkg.test()` forces `--check-bounds=yes` (see `test_threading.jl:25-27`).
This is deliberate, not overhead to shave off: the code indexes ghost-cell
buffers by hand under `@inbounds`/`@turbo`, exactly where bounds checking
earns its cost. Do not disable or "optimize" it away.

Shared fixtures live in `test/testutils.jl` (`synthetic_fields`, `quiet`,
`with_tempdir`, `gmean`). Prefer exact `==` over `isapprox` where the claim
is bit-identity — see the `-t 1` vs `-t 2` monthly-mean comparison in
`test_threading.jl`. Do not run the model to test something the model does
not do.
