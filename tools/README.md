# Maintainer tools

Scripts for maintainers. None of them is needed to use the package. Run them
from the repository root with `julia --project=. tools/<folder>/<script>.jl`.

| Folder | Purpose | Script | What it does |
|:-------|:--------|:-------|:-------------|
| `dataset/` | Build and publish the model's input dataset | `fetch_era5_data.py` | Fetches an ERA5 climatology and writes GREB `.bin` files |
| | | `convert_greb_to_jld2.jl` | Converts GREB `.bin` files into `greb_input_data/` |
| | | `package_dataset.jl` | Builds the dataset archive and its SHA256 for the DataDep |
| `validation/` | Check the model against a reference | `bit_identity.jl` | Saves every record field of fixed runs, then compares a later build with exact equality; for refactors that must not change results |

Where a new script goes:

| If it... | Folder |
|:---------|:-------|
| produces or publishes `greb_input_data/` | `dataset/` |
| turns an external record into a forcing-series file the model reads | `forcing/` |
| runs the model to measure a property of the model | `diagnostics/` |
| runs the model and compares with observed data | `validation/` |

Scripts that run the model should read the local dataset with
`greb_data_dir(; allow_download=false)` and never download it.

External data (NetCDF sources, converted records, observations) lives in the
gitignored `Data/` folder. Converters that need a package the model does not
depend on, such as NCDatasets, say so in their header; add it to your default
environment rather than to the project.
