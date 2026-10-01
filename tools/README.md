# Maintainer tools

Scripts for maintainers. None of them is needed to use the package. Run them
from the repository root with `julia --project=. tools/<folder>/<script>.jl`.

| Folder | Purpose | Script | What it does |
|:-------|:--------|:-------|:-------------|
| `dataset/` | Build and publish the model's input dataset | `fetch_era5_data.py` | Fetches an ERA5 climatology and writes GREB `.bin` files |
| | | `convert_greb_to_jld2.jl` | Converts GREB `.bin` files into `greb_input_data/` |
| | | `package_dataset.jl` | Builds the dataset archive and its SHA256 for the DataDep |
| `diagnostics/` | Measure a property of the model | `budget.jl` | Runs a preset with a per-step observer and prints whether each store changes by the sum of its flows, how often a limiter held a cell, and the global-mean flows |
| `validation/` | Check the model against a reference | `bit_identity.jl`, `preset_reference.jl` | `bit_identity.jl` saves every record field of every experiment preset, then compares a later build with exact equality; for refactors that must not change results. `preset_reference.jl` writes `test/data/preset_reference.jl`, what each preset imposes, which `test/test_presets.jl` checks |

Where a new script goes:

| If it... | Folder |
|:---------|:-------|
| produces or publishes `greb_input_data/` | `dataset/` |
| runs the model to measure a property of the model | `diagnostics/` |
| runs the model and compares with observed data | `validation/` |

Scripts that run the model should read the local dataset with
`greb_data_dir(; allow_download=false)` and never download it.

External data (NetCDF sources, converted records, observations) lives in the
gitignored `Data/` folder. Converters that need a package the model does not
depend on, such as NCDatasets, say so in their header; add it to your default
environment rather than to the project.
