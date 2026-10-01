---
paths: ["src/**/*.jl"]
---

`circulation!` is ~93-95% of per-timestep cost after the ghost-cell change
(`865ae01`), down from ~98% before it (re-measured 2026-09-22: 94.5%) —
anything here outweighs everything else in `src/` combined.

Longitude buffers carry `nghost = 3` periodic ghost cells at each end
(`xghost = xdim + 2*nghost`, `core/constants.jl`); `nghost` equals the
stencil reach, so widening one means widening both.

`CircularArrays.jl` is not an option here: `LoopVectorization.check_args`
rejects it and `@turbo` silently falls back to a scalar loop rather than
erroring — confirmed via `code_native`, ~12 `gather`-family instructions in
the working ghosted version's baseline vs. a fully scalar inner loop
(`vmovsd`/`vmulsd`/...) for the circular version. A ghosted `Matrix{Float32}`
does not have this failure mode. When touching a hot loop here, check
`code_native` for `gather` (should be 0 in `_diffusion!`/`_advection!`) and
for `cvtss2sd`/`cvtsd2ss` (should be 0 — a stray `Float64` literal).

**Bit-identical is the bar** for changes to `circulation!`. One accepted
deviation exists: since the ghost-cell change `circulation!(q)` differs in 56
of 4608 cells (max 7.3e-7 relative) because `@turbo` contracts FMAs
differently on packed loads; do not add a second one without the same rigor.
Optimization wins do not add up — measure the combination, not each piece in
isolation. Use generic indexing (`axes`, `eachindex`) rather than
assuming `1:n`; `@inbounds` is deliberate throughout this file and is not to
be added — or removed — casually.
