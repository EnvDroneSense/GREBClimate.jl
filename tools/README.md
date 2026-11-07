# Maintainer tools

Scripts for maintainers. None of them is needed to use the package. Run them
from the repository root with `julia --project=. tools/<folder>/<script>.jl`
(the one Python script with `python`).

| Folder | Purpose | Script | What it does | Page |
|:-------|:--------|:-------|:-------------|:-----|
| `dataset/` | Build and publish the model's input dataset | `fetch_era5_data.py` | Fetches an ERA5 climatology and writes GREB `.bin` files | [Dataset tools](../docs/src/dev/dataset_tools.md) |
| | | `convert_greb_to_jld2.jl` | Converts GREB `.bin` files and the CO2 text tables into `greb_input_data/` | |
| | | `package_dataset.jl` | Builds one `tar.xz` per pack and prints the `DATA_PACKS` block | |
| | | `fields.toml`, `fields.jl` | Not run: the file list, pack and metadata of each field, and the reader of that list | |
| `diagnostics/` | Measure a property of the model | `budget.jl` | Reports whether each store closes, how often a limiter held a cell, and the global-mean flows | [Diagnostics](../docs/src/dev/diagnostics.md) |
| `experiments/` | Run the model to answer a question | `sensitivity.jl`, `co2_sweep.jl` | The effective climate sensitivity of an abrupt-CO2 run, and the last-year climate over a range of CO2 levels | [Experiments](../docs/src/dev/experiments.md) |
| | | `helpers.jl` | The module `ExperimentTools`: annual means, the Gregory estimate and map diagnostics | |
| `validation/` | Check the model against a reference | `bit_identity.jl` | Saves every record field of every preset, then compares a later build with exact equality | [Validation](../docs/src/dev/validation.md) |

Building and publishing the dataset (details in the `publish-dataset` skill):

```bash
julia --project=. tools/dataset/convert_greb_to_jld2.jl Data greb_input_data
julia --project=. tools/dataset/package_dataset.jl greb_input_data dataset_packs
```

The packs are the unit of release: edit `DATA_PACKS` in `src/data.jl` with the printed block,
upload all five archives to the release tag in `DATA_RELEASE_TAG`, and raise `DATASET_VERSION`
(in `src/data.jl` and in `fields.jl`) whenever the format of the files changes. Tests of the
scripts, and of the agreement between the file list, the packs and the loaders, are in
`test/test_dataset_tools.jl`. They are not part of the full suite; run them by name with
`Pkg.test(test_args=["dataset_tools"])` after changing the tools, `fields.toml` or the loaders.

Where a new script goes:

| If it... | Folder |
|:---------|:-------|
| produces or publishes `greb_input_data/` | `dataset/` |
| runs the model to measure a property of the model | `diagnostics/` |
| runs the model many times, or reduces its output, to answer a question | `experiments/` |
| compares the model with a reference: observations, a saved snapshot or a stored table | `validation/` |

Scripts that run the model should read the local dataset with
`greb_data_dir(; allow_download=false)` and never download it.

External data (NetCDF sources, converted records, observations) lives in the
gitignored `Data/` folder. Converters that need a package the model does not
depend on, such as NCDatasets, say so in their header; add it to your default
environment rather than to the project.
