---
name: docs-check
description: Verify GREBClimate.jl's documentation matches the code - run every Julia snippet in README.md and docs/src/ against the real package, and check file references, exported names and version claims. Use before claiming docs are updated, after renaming or changing any exported function's signature, and when reviewing changes to README.md, DATA_README.md or docs/.
---

# Docs check

The README's quick-start once ran the model on an all-zero climatology and reported a -40 degC world as success, because the snippet discarded the return value of `load_climatology`. It exited 0. Reading a snippet is not checking it: **run it.**

## When to use this

- Before telling the user the docs are correct or updated.
- After changing the signature, name or return value of anything exported from `GREBClimate`.
- When reviewing edits to `README.md`, `DATA_README.md`, `CONTRIBUTING.md` or `docs/src/`.
- After a rename. A stale `docs/Manifest.toml` left by the `GREB` to `GREBClimate` rename broke the docs build for weeks.

## The checks

| # | Check | Command |
|---|---|---|
| 1 | Run the documented snippets | See below |
| 2 | Local Markdown links resolve, relative to the linking file | `julia .claude/skills/docs-check/check_docs.jl` |
| 3 | Julia version claims agree across `README.md`, `docs/src/index.md` and the `Project.toml` compat bound | same |
| 4 | Names in `[`name`](@ref)` are exported | same |
| 5 | No dead `.jld2` file (a NOTE, not a failure; the check expands `$suffix` in the ENSO file names) | same |
| 6 | Documenter cross-references, docstrings of exports, doctests | `julia --project=docs docs/make.jl` |

The static script uses only the standard library, does not load the package, runs in seconds and exits nonzero on failure. The files it checks are listed in `DOC_FILES`; add a new page there. Do not replace it with an ad-hoc `grep`: that resolves every path from the repository root and gives about 70 false positives.

Check 6 is set to `checkdocs = :exports`, so an exported name without a docstring fails. A clean build prints only "could not auto-detect the building environment. Skipping deployment." and a size warning for `api.md`. If `Pkg.develop` says the UUID already exists, delete the gitignored `docs/Manifest.toml`; it regenerates.

## Running the snippets

The same flow is documented in `README.md` (quick-start and "Running the model"), `docs/src/tutorial.md`, `docs/src/examples.md` and `examples/quickstart.jl`, and they drift independently.

```bash
julia --project=. -t 2,0 examples/quickstart.jl      # the canonical path
julia --project=. examples/quickstart.jl --data-free  # no dataset; numbers mean nothing
```

For a snippet the example does not cover, paste it verbatim into a scratch file and run it. Verbatim matters: the bug above reproduces only if the missing `fields=` stays missing.

Check the numbers, not the exit code. A `:full_model` control run on the real dataset with the default `SpinUp(3)`:

| | Healthy | On a zero climatology |
|---|---|---|
| The logged line | `1970: Ts global mean 13.82 degC` | `Ts global mean NaN degC` |
| `mean(result.ctrl[1].Ts)` | 276.94 K | NaN |

Measured 2026-09-29 and unchanged since. A line above 14 degC means the spin-up did not run: the control then used corrections that do not fit its configuration and drifted warm. The logged value is the area-weighted annual mean and the array mean is an unweighted mean of one month, so compare each only with its own column.

## Reporting

Say which checks ran and what they printed. "Docs look right" is not a result; "ran the quick-start, global-mean Ts 13.82 degC, docs build clean" is. If a check fails, fix the docs rather than the check, unless the code is what is wrong.
