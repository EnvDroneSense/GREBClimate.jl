# Maintainer tools

Scripts for maintainers. None of them is needed to use the package. Run them
from the repository root with `julia --project=. tools/<folder>/<script>.jl`
(the one Python script with `python`).

| Folder | Purpose | Script | What it does | Page |
|:-------|:--------|:-------|:-------------|:-----|
| `dataset/` | Build and publish the model's input dataset | `fetch_era5_data.py` | Experimental, not in use. Fetches an ERA5 climatology and writes GREB `.bin` files | [Dataset tools](../docs/src/dev/dataset_tools.md) |
| | | `convert_greb_to_jld2.jl` | Converts GREB `.bin` files and the CO2 text tables into `greb_input_data/` | |
| | | `package_dataset.jl` | Builds one `tar.xz` per pack and prints the `DATA_PACKS` block | |
| | | `fields.toml`, `fields.jl` | Not run: the file list, pack and metadata of each field, and the reader of that list | |
| `checks/` | Measure the model, or check it against a reference | `budget.jl` | Reports whether each store closes, how often a limiter held a cell, and the global-mean flows | [Diagnostics](../docs/src/dev/diagnostics.md) |
| | | `bit_identity.jl` | Saves every record field of every preset, then compares a later build with exact equality | [Validation](../docs/src/dev/validation.md) |

Two maintainer scripts live next to what they serve, not here:

| Script | What it does | Page |
|:-------|:-------------|:-----|
| `test/support/preset_reference.jl` | Regenerates `test/support/data/preset_reference.jl` after a preset is changed on purpose | [Testing](../docs/src/dev/testing.md) |
| `.claude/skills/docs-check/check_docs.jl` | Static checks of the Markdown docs: local links and Julia version claims | [Code conventions](../docs/src/dev/conventions.md) |

The experiment scripts (climate sensitivity, CO2 sweep) and their helpers are not here. They are in the separate repository [GREBExperiments](https://github.com/EnvDroneSense/GREBExperiments).

Building and publishing the dataset (details in the `publish-dataset` skill):

```bash
julia --project=. tools/dataset/convert_greb_to_jld2.jl Data greb_input_data
julia --project=. tools/dataset/package_dataset.jl greb_input_data dataset_packs
```

The packs are the unit of release: edit `DATA_PACKS` in `src/data.jl` with the printed block,
upload all five archives to the release tag in `DATA_RELEASE_TAG`, and raise `DATASET_VERSION`
(in `src/data.jl` and in `fields.jl`) whenever the format of the files changes. The scripts have
no tests; run them by hand at a release. `test/test_io.jl` checks that the file list, the packs
and the version agree with the loaders.

Where a new script goes:

| If it... | Folder |
|:---------|:-------|
| produces or publishes `greb_input_data/` | `dataset/` |
| runs the model to measure a property, or compares it with a reference such as a saved snapshot | `checks/` |

`ToolBox.jl` holds what several scripts share: the local dataset lookup, a quiet run, and observers that keep
one result per phase or sum area-weighted terms. Load it with
`include(joinpath(@__DIR__, "..", "ToolBox.jl")); using .ToolBox`.

Scripts that run the model should read the local dataset with
`greb_data_dir(; allow_download=false)` and never download it.

External data (NetCDF sources, converted records, observations) is in the
`Data/` folder, which Git ignores. Converters that need a package the model does not
depend on, such as NCDatasets, say so in their header; add it to your default
environment rather than to the project.
