# The dataset's file list, its packs and its metadata, shared by the converter and
# the packager. The per-field data (description, units, source, period, pack) is in
# fields.toml. The list is the tools' own, so a file can be added there before the
# model reads it; test/test_io.jl checks that every file the loaders read is in it,
# in the pack the loaders name for it.

using TOML

# Format version written into every file; GREBClimate.DATASET_VERSION is the
# version the loaders require (a test checks they are equal).
const DATASET_VERSION = 2

# Files that hold several fields; the converter writes them from special-case code.
const FLUX_CORRECTION_NAMES = ("Tsurf_flux_correction", "vapour_flux_correction", "Tocean_flux_correction")
const COMBINED_FILE_NAMES = ("flux_corrections", "ipcc_scenarios", "solar_paleo", "solar_eccentricity",
    "solar_obliquity", "historical_emissions_population")

const FIELDS = TOML.parsefile(joinpath(@__DIR__, "fields.toml"))

# Every single-field file the model reads
const MODEL_FIELD_NAMES = Set{String}(keys(FIELDS["field"]))

# ---- packs: `core` is what a default run needs, the others are optional downloads ----

const PACK_NAMES = (:core, :era, :cmip5, :enso, :solar)

"The pack that holds the file at `relpath` of the dataset tree (`/` as separator); combined files are `core`."
function pack_of(relpath::AbstractString)
    startswith(relpath, "solar_scenarios/") && return :solar
    entry = get(FIELDS["field"], first(splitext(basename(relpath))), nothing)
    return entry === nothing ? :core : Symbol(entry["pack"])
end

# ---- metadata ----

# The grid of every 2D and 3D field: 96 x 48 points of 3.75 degrees, and a mean
# year of 730 steps of 12 h (time in days from 1 January 00:00).
const GRID_COORDS = (lon=collect(range(0.0, step=3.75, length=96)),
    lat=collect(range(-88.125, step=3.75, length=48)),
    time=collect(range(0.0, step=0.5, length=730)))

const SOURCE_GREB = FIELDS["sources"]["greb"]
const SOURCE_IPCC = FIELDS["sources"]["ipcc"]

"The `(description, units, source, period)` of the field file `name` (without ending)."
function field_info(name::AbstractString)
    entry = get(FIELDS["field"], name, nothing)
    entry === nothing && error("no metadata for field $name: add it to tools/dataset/fields.toml")
    return (description=entry["description"], units=entry["units"], source=FIELDS["sources"][entry["source"]],
        period=entry["period"])
end

# The same for the three flux corrections, by name
const FLUX_CORRECTION_INFO = NamedTuple(Symbol(k) => (description=FIELDS["flux_correction"][k]["description"],
                                            units=FIELDS["flux_correction"][k]["units"])
for k in FLUX_CORRECTION_NAMES)
