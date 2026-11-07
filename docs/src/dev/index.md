# Developer tools

The scripts in `tools/` are for maintainers. Nobody needs them to use the package. Run each from the repository root with `julia --project=. tools/<folder>/<script>.jl`. This section has one page per group of tools; a group gets a page here once it needs more explanation than its table row in `tools/README.md`.

| Folder | Purpose | Page |
|:-------|:--------|:-----|
| `tools/dataset/` | Build and publish the input dataset | [Dataset tools](dataset_tools.md) |
| `tools/diagnostics/` | Measure a property of the model, such as its budgets | `tools/README.md` |
| `tools/experiments/` | Run the model to answer a question, such as its climate sensitivity | `tools/README.md` |
| `tools/validation/` | Compare the model with a reference | `tools/README.md` |

Scripts that run the model read the local dataset with `greb_data_dir(; allow_download=false)` and never download it.
