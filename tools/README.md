# Maintainer tools

Scripts for maintainers. None of them is needed to use the package. Run them
from the repository root with `julia --project=. tools/<folder>/<script>.jl`.
The scripts that run the model read the local dataset through `common.jl`, and
never download it.

| Folder | Purpose | Script | What it does |
|:-------|:--------|:-------|:-------------|
| `dataset/` | Build and publish the model's input dataset | `fetch_era5_data.py` | Fetches an ERA5 climatology and writes GREB `.bin` files |
| | | `convert_greb_to_jld2.jl` | Converts GREB `.bin` files into `greb_input_data/` |
| | | `package_dataset.jl` | Builds the dataset archive and its SHA256 for the DataDep |
| `forcing/` | Convert published forcing records to the forcing-series format | `convert_aerosol_to_jld2.jl` | GISS stratospheric-aerosol NetCDF to an `aerosol_aod` file |
| `diagnostics/` | Measure how the model responds to a forcing | `aerosol_forcing_per_aod.jl` | Shortwave change per unit optical depth; the source of `tau_scale` |
| | | `shortwave_pulse_response.jl` | Temperature response to a one-year and a step shortwave forcing |
| | | `aerosol_profile_check.jl` | Kernel optical depth against the aerosol records and Stothers (1996), per latitude band, for reference eruptions |
| `validation/` | Compare the model with observations | `validate_volcanic.jl` | Agung, El Chichón and Pinatubo against GISTEMP and MEI |

Where a new script goes:

| If it... | Folder |
|:---------|:-------|
| produces or publishes `greb_input_data/` | `dataset/` |
| turns an external record into a forcing-series file the model reads | `forcing/` |
| runs the model to measure a property of the model | `diagnostics/` |
| runs the model and compares with observed data | `validation/` |

External data (NetCDF sources, converted records, observations) lives in the
gitignored `Data/` folder. Converters that need a package the model does not
depend on, such as NCDatasets, say so in their header; add it to your default
environment rather than to the project.
