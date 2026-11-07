# Dataset tools

These tools turn the raw GREB files into the `.jld2` dataset and package it as release assets. Users never run them: they get the finished packs through `greb_data_dir()` (see [Input data](../data.md)).

```
Data/ (raw .bin, .ctl, .txt)  --convert-->  greb_input_data/ (39 .jld2)  --package-->  dataset_packs/ (5 tar.xz)  --upload-->  GitHub release data-v2
```

| File | Role |
|:-----|:-----|
| `tools/dataset/fields.toml` | One line per field, with its description, units, source, period and pack. Add a new field here. |
| `tools/dataset/fields.jl` | Reads `fields.toml` and answers which files exist and which pack each belongs to. |
| `tools/dataset/convert_greb_to_jld2.jl` | Converts the raw files to `.jld2`. Each file gets its metadata and `dataset_version`. |
| `tools/dataset/package_dataset.jl` | Checks the converted tree, builds one `tar.xz` per pack and prints the `DATA_PACKS` block. |
| `tools/dataset/fetch_era5_data.py` | Work in progress. It fetches ERA5 fields, but nothing in the dataset uses its output yet. |

## Packs

| Pack | Files | Archive | Holds | Needed by |
|:-----|------:|--------:|:------|:----------|
| `core` | 17 | 103 MB | static fields, NCEP climatology, common fields, flux corrections, solar radiation, CO2 tables | every run |
| `era` | 4 | 38 MB | ERA-Interim climatology | `dataset = :era` |
| `cmip5` | 5 | 54 MB | CMIP5 RCP8.5 anomalies | the `:rcp85` experiment |
| `enso` | 10 | 121 MB | El Nino and La Nina anomalies | `:elnino`, `:lanina` |
| `solar` | 3 | 4 MB | paleo and orbital solar tables | `:paleo_231kyr`, `:eccentricity`, `:obliquity` |

A file belongs to the pack named on its line in `fields.toml`. When a file is missing, the loaders take the pack to suggest from `_OPTIONAL_PACK_FILES` in `src/data.jl`. The dataset tools test checks that the two lists agree.

## Build

```bash
julia --project=. tools/dataset/convert_greb_to_jld2.jl Data greb_input_data
julia --project=. tools/dataset/package_dataset.jl greb_input_data dataset_packs
```

The converter stops at the first failure. It skips raw files that the model never reads; `--all` converts them for inspection, but the packager rejects such a tree. The packager checks the tree before it builds anything. Every file the model reads must be present, no other file may be, and each file must carry its metadata. A correct build prints `39 .jld2 files, 438.8 MB unpacked` and the pack sizes above. Any other count means the tree is wrong.

The archives are reproducible. Entries are sorted, owner and group are zeroed, the timestamp is pinned and `xz` runs on one thread, so the same input always gives the same checksums. Never change the pinned `--mtime`. If a plain regeneration changes a checksum, either the data changed or the GNU `tar` or `xz` version differs, so find out which before you use the new value. The checksum that counts is the one of the file you upload.

## Publish

Three values must agree. Otherwise users get a 404 or a checksum failure that appears only on a machine without a local copy.

| Value | Where | Must match |
|:------|:------|:-----------|
| release tag | `DATA_RELEASE_TAG` in `src/data.jl` | the git release tag |
| asset names | `DATA_PACKS[pack].archive` | the uploaded file names, exactly |
| hashes and sizes | `DATA_PACKS[pack].sha256` and `.bytes` | the built archives |

Paste the printed `DATA_PACKS` block over the old one, then upload. `gh` takes the asset name from the local file name, so upload the files under the names the packager gave them.

```bash
gh release create data-v2 dataset_packs/greb_input_data-v2-*.tar.xz --repo EnvDroneSense/GREBClimate.jl \
  --title "GREB input dataset v2" --notes "JLD2 input dataset, version 2: five packs."
# to replace the assets of an existing release:
# gh release upload data-v2 dataset_packs/*.tar.xz --clobber --repo EnvDroneSense/GREBClimate.jl
```

When the file format changes, raise `DATASET_VERSION` in `src/data.jl` and in `tools/dataset/fields.jl`, and use a new tag such as `data-v3`. The DataDep name contains the version, so an older cached dataset is never mistaken for the new one. Do not replace assets that older installs may have cached.

## Verify

The dataset tools test does not run with the full suite. Run it by name after you change the tools, `fields.toml` or the loaders in `src/io.jl` and `src/data.jl`.

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["io", "dataset_tools"])'
for p in core era cmip5 enso solar; do
  curl -sIL -o /dev/null -w "$p %{http_code}\n" "$(julia --project=. -e "using GREBClimate; print(GREBClimate._pack_url(:$p))")"
done
```

A `200` for each pack means the assets are reachable. The download path only runs when no local dataset exists. Test it in a clean checkout, a container or a temporary depot (`JULIA_DEPOT_PATH`), with `GREB_DATA` unset and `DATADEPS_ALWAYS_ACCEPT=true`, by calling `greb_data_dir()`. Do not move `greb_input_data/` aside to simulate this: a rename-and-restore once destroyed the directory.

## Failures

| Symptom | Cause |
|:--------|:------|
| DataDeps 404 | The release or asset is missing, or the asset name differs from `DATA_PACKS[pack].archive`. |
| DataDeps hash mismatch | `DATA_PACKS` was not updated, or the asset was rebuilt after its hash was recorded. |
| Process hangs without output | A consent prompt is waiting on stdin. Set `DATADEPS_ALWAYS_ACCEPT=true`. |
| Packager reports extra files | The tree was built with `--all`. |
| Packager reports missing files | The conversion was incomplete, or `fields.toml` lists a field the converter cannot produce. |
| Packager reports metadata problems | The tree came from an older converter. Regenerate it. |
| `ArgumentError` about dataset version 1 | The dataset predates the metadata release. Download the current one. |
| `ArgumentError` about an optional pack | The file is in an optional pack. Call `download_greb_pack(:era)` or the pack it names. |
