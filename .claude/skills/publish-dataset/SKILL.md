---
name: publish-dataset
description: Build and publish the GREB JLD2 input dataset as GitHub release assets (one tar.xz per pack), and reconcile DATA_PACKS in src/data.jl. Use after regenerating greb_input_data/ from raw .bin files, when the DataDep download 404s or fails checksum verification, or when asked to cut/update a dataset release.
---

# Publish the dataset

The full reference (packs, reproducibility, failure table) is `docs/src/dev/dataset_tools.md`. This is the checklist.

1. Convert and build. Expect `39 .jld2 files, 438.8 MB unpacked` and packs of 17, 4, 5, 10 and 3 files:

   ```bash
   julia --project=. tools/dataset/convert_greb_to_jld2.jl Data greb_input_data
   julia --project=. tools/dataset/package_dataset.jl greb_input_data dataset_packs
   ```

2. Paste the printed `DATA_PACKS` block over the old one in `src/data.jl`. A changed checksum after a plain regeneration is not expected: investigate before pasting. Never touch the `--mtime` pin.
3. Upload the archives, named exactly as `DATA_PACKS[pack].archive`, to the release `DATA_RELEASE_TAG`:

   ```bash
   gh release create data-v2 dataset_packs/greb_input_data-v2-*.tar.xz --repo EnvDroneSense/GREBClimate.jl \
     --title "GREB input dataset v2" --notes "JLD2 input dataset, version 2: five packs."
   ```

4. Verify: `Pkg.test(test_args=["io", "dataset_tools"])` (the dataset tools test is not in the full suite; run it by name), a `curl -sIL` that returns 200 for each `_pack_url(pack)`, and a download on a machine without a local dataset (`GREB_DATA` unset, `DATADEPS_ALWAYS_ACCEPT=true`).

An incompatible format change needs a new tag (`data-v3`) and a raised `DATASET_VERSION` in `src/data.jl` and `tools/dataset/fields.jl`. Never replace assets older installs may have cached, and never move `greb_input_data/` aside to test the download.
