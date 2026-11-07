# Validation

`tools/validation/bit_identity.jl` answers one question: did a change alter any number the model returns? It is for refactors and optimisations that claim "no change to results". It needs the local dataset.

```bash
julia --project=. -t 2,0 tools/validation/bit_identity.jl save    <snapshot.jld2>
julia --project=. -t 2,0 tools/validation/bit_identity.jl compare <snapshot.jld2>
```

`save` runs every experiment preset plus two cases with switches off (44 cases) for one year of spin-up, one control year and one scenario year, and stores every monthly record field. `compare` runs the same cases again and compares with exact equality, where `NaN` equals `NaN`. A snapshot is about 290 MB and a run takes a few minutes.

```text
Threads.nthreads() = 2, Julia 1.12.6
  ran a1b
  ran co2_10x
  ...
  ran sst_plus1
IDENTICAL: 0 differing values in 1320 arrays
```

| Output | Meaning |
|:-------|:--------|
| `IDENTICAL: 0 differing values` | No stored value changed. A refactor may claim "no change to results" |
| `DIFFERENT: N differing values` | One line per array before it: `<case>/<phase>/<field>: n differing values, max absolute difference`. Start with the first case. The exit code is 1 |
| `new cases` or `new record fields (not compared)` | A preset or field the snapshot lacks. Not a difference; take a new snapshot to include it |
| A warning about the Julia version | The snapshot is from another Julia release, where exact equality is not expected |
| An error about threads | The snapshot was taken with another thread count. Rerun with the `-t` it names |
| An error about missing arrays | A preset was renamed or removed. Take the snapshot again on the commit before the change |

## Workflow

| Step | Command |
|:-----|:--------|
| 1 | On the commit before the change, run `save` with `-t 1` |
| 2 | On the commit after the change, run `compare` with `-t 1` |
| 3 | Repeat both with `-t 2,0` and a second snapshot |

The thread count matters because the two circulation calls run concurrently only above one thread. Snapshots compare only on the same machine, Julia version and thread count, because `@turbo` gives different results on different CPUs and releases. They are not committed.

## How it differs from other checks

| Check | Compares | Tolerance |
|:------|:---------|:----------|
| `bit_identity.jl` | Your change against the commit before it, on every preset | None |
| `test/test_golden.jl` | The model against stored reference numbers | A tolerance, and it skips in CI |
| [Fortran comparison](fortran.md) | The model against the original Fortran | The limits of the gate |
| [Diagnostics](diagnostics.md) | The model's budgets with themselves | Rounding |

A change that is meant to alter results does not pass `bit_identity.jl` by design. For that, check the golden test and say in the pull request which numbers moved.
