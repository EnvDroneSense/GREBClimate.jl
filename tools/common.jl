# Shared by the Julia scripts in tools/ that run the model. Include after
# `using GREBClimate`.

const REPO = normpath(joinpath(@__DIR__, ".."))

# The local dataset; never downloads.
const DATA_DIR = something(greb_data_dir(; allow_download=false),
                           joinpath(REPO, "greb_input_data"))
isdir(DATA_DIR) || error("dataset not found at $DATA_DIR")
