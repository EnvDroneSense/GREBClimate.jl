# =============================================================================
# data.jl - locating the JLD2 input dataset (not shipped with the package).
# `greb_data_dir()` resolves where it lives, falling back to a DataDeps.jl
# download of the core pack; `download_greb_pack` adds an optional pack.
# =============================================================================

"""
Version of the dataset format this package reads. Every dataset file carries
its `dataset_version`; the loaders refuse a dataset of another version.
"""
const DATASET_VERSION = 2

"""
Name of the DataDep of the core pack, and the folder name it is cached under
inside DataDeps' storage (`~/.julia/scratchspaces/<uuid>/datadeps/` by default).
The version is part of the name, so a dataset downloaded for an older format
is not mistaken for this one.
"""
const DATA_DEP_NAME = "GREB-input-data-v$DATASET_VERSION"

"""
Release tag, and for each pack the asset name, size and SHA256 of the published
archive. The packs are built from a validated `greb_input_data/` tree with
`tools/dataset/package_dataset.jl`, which also prints the `DATA_PACKS` block.

`core` is what a default run needs. The optional packs add: `era` the
ERA-Interim climatology (`dataset = :era`), `cmip5` the CMIP5 RCP8.5 anomalies,
`enso` the El Nino and La Nina anomalies, `solar` the paleo and orbital solar
tables.
"""
const DATA_RELEASE_TAG = "data-v2"
const DATA_PACKS = (
    core=(archive="greb_input_data-v2-core.tar.xz", bytes=102955736,
        sha256="fc3892f411945117dc08ff3e79bd1b99e0e2a5d67ee4dad1994d36c137a84147"),
    era=(archive="greb_input_data-v2-era.tar.xz", bytes=38378900,
        sha256="043c46a0964b6747ad01d18743681c9f5143de7a5198b3ed79fb5320c865f8f6"),
    cmip5=(archive="greb_input_data-v2-cmip5.tar.xz", bytes=54359544,
        sha256="f36e9fa819394e1748587782d1158f6ea1024c30443ed5ca3bf10335d4c75491"),
    enso=(archive="greb_input_data-v2-enso.tar.xz", bytes=120905740,
        sha256="e5ff850ed0672de51bbbbafccb62ae0b5fac400a5d1c0fbde9d14665b13f61a4"),
    solar=(archive="greb_input_data-v2-solar.tar.xz", bytes=4440816,
        sha256="9bf3593a9af0223ee8d73e4f72674b4439243b046d18df741d14e0d7c642253e")
)

const _OPTIONAL_PACKS = (:era, :cmip5, :enso, :solar)

_pack_dep_name(pack::Symbol) = pack === :core ? DATA_DEP_NAME : "$DATA_DEP_NAME-$pack"
_pack_url(pack::Symbol) = "https://github.com/EnvDroneSense/GREBClimate.jl/releases/download/" *
                          "$DATA_RELEASE_TAG/$(DATA_PACKS[pack].archive)"

const _PACK_CONTENT = (
    core="the static fields, the NCEP climatology, the fields common to both datasets, the flux corrections, " *
         "the solar radiation and the CO2 tables",
    era="the ERA-Interim climatology (dataset = :era)",
    cmip5="the CMIP5 RCP8.5 anomalies (the :rcp85 experiment)",
    enso="the El Nino and La Nina anomalies (the :elnino and :lanina experiments)",
    solar="the paleo and orbital solar tables (the :paleo_231kyr, :eccentricity and :obliquity experiments)"
)

function _pack_message(pack::Symbol)
    megabytes = round(Int, DATA_PACKS[pack].bytes / 1e6)
    optional = round(Int, sum(DATA_PACKS[p].bytes for p in _OPTIONAL_PACKS) / 1e6)
    follow = pack === :core ? "
The four optional packs ($optional MB) are not part of this download. A run that needs one stops with a message, and download_greb_pack(pack) then asks for it.
" : ""
    return """
    The GREB input dataset, pack "$pack" (~$megabytes MB download): $(_PACK_CONTENT[pack]).
    $follow

    Derived from the input files of the original GREB model by Dietmar Dommenget
    and colleagues (Monash University), converted to JLD2. It combines fields from
    several upstream sources:

      NCEP/NCAR reanalysis      surface temperature, winds, humidity, soil moisture
      ERA-Interim (ECMWF)       surface temperature, winds, humidity, omega
      ISCCP                     cloud cover
      WOCE                      ocean mixed-layer depth
      CMIP5 (RCP8.5 ens. mean)  climate-change anomaly forcing
      IPCC                      RCP/SSP CO2 scenario tables

    """
end

# The files of each optional pack, as paths in the dataset directory. The dataset
# tools place a file in a pack by its name (`pack_of` in tools/dataset/fields.jl);
# test/test_io.jl checks that the two agree.
const _OPTIONAL_PACK_FILES = (
    era=["climatology/$name.jld2" for name in values(_CLIMATOLOGY_FILES.era) if !startswith(name, "ncep.")],
    cmip5=["climatology/$name.jld2" for name in values(_boundary_anomaly_files(:cmip5_rcp85))],
    enso=["climatology/$name.jld2" for source in (:elnino, :lanina) for name in values(_boundary_anomaly_files(source))],
    solar=["solar_scenarios/solar_$kind.jld2" for kind in (:paleo, :eccentricity, :obliquity)]
)

# The optional pack that holds the file `name` (with or without ".jld2"), or `nothing`
function _pack_of_file(name::AbstractString)
    file = endswith(name, ".jld2") ? name : name * ".jld2"
    for pack in _OPTIONAL_PACKS
        any(f -> basename(f) == file, getfield(_OPTIONAL_PACK_FILES, pack)) && return pack
    end
    return nothing
end

"""
    greb_data_dir(path = nothing; allow_download = true, packs = :all) -> String or nothing

Return the directory holding the JLD2 input dataset, resolving in this order:

1. `path`, if given and non-empty.
2. `ENV["GREB_DATA"]`, if set.
3. `greb_input_data/` in the package's top folder, if it exists.
4. An already-downloaded `$DATA_DEP_NAME` DataDep, if its cache is present.
5. The `$DATA_DEP_NAME` DataDep - downloading it on first use, after asking.

Pass `allow_download = false` to stop after step 4 and return `nothing` when no
dataset is available locally.

The first download (step 5) asks once: "Download the whole dataset (321 MB)?".
Yes fetches all five packs without further questions. No fetches the core pack
only (103 MB), which is enough for a default run, after DataDeps' usual question
for that pack; a run that later needs an optional pack (the ERA climatology, the
anomalies of an experiment, the orbital solar tables) stops with an error that
names the pack, and [`download_greb_pack`](@ref) asks for and adds it. With
`ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"` nothing is asked and everything is
fetched. `packs = :core` skips the question and fetches the core pack only.
A dataset that is already downloaded, or in a folder of your own (steps 1 to 4),
is never extended.

The result is a plain path, suitable for [`load_climatology`](@ref) and
`greb_model!`'s `jld2_dir`:

```julia
dir    = greb_data_dir()
fields = load_climatology(dir; dataset = :ncep)
result = greb_model!(RunSpec(), preset(:co2_double); jld2_dir = dir, fields = fields)
```
"""
function greb_data_dir(path::Union{Nothing,AbstractString}=nothing;
        allow_download::Bool=true, packs::Symbol=:all)
    packs in (:all, :core) ||
        throw(ArgumentError("packs must be :all or :core, got :$packs"))
    if path !== nothing && !isempty(path)
        isdir(path) || error("GREB dataset directory not found: $path")
        return String(path)
    end

    env = get(ENV, "GREB_DATA", "")
    if !isempty(env)
        isdir(env) || error("GREB_DATA points at a directory that does not exist: $env")
        return env
    end

    local_dir = normpath(joinpath(@__DIR__, "..", "greb_input_data"))
    isdir(local_dir) && return local_dir

    # A dataset downloaded earlier is used even with `allow_download = false`:
    # reading it needs no network.
    cached = _cached_datadep_path()
    cached === nothing || return cached

    allow_download || return nothing

    # Resolved lazily: this is the only branch that can trigger a download, and
    # it must never run at precompile time
    if packs === :all && _download_everything()
        # One yes covers the five packs: DataDeps is told to accept them all
        return withenv("DATADEPS_ALWAYS_ACCEPT" => "true") do
            dir = @datadep_str DATA_DEP_NAME
            download_greb_pack(:all; dir)
        end
    end
    return @datadep_str DATA_DEP_NAME
end

# The one question of a first download: everything (true), or the core pack only.
# Without an answer (no stdin) it is the core pack, whose own prompt then follows.
function _download_everything()
    lowercase(get(ENV, "DATADEPS_ALWAYS_ACCEPT", "")) == "true" && return true
    total = round(Int, sum(p.bytes for p in values(DATA_PACKS)) / 1e6)
    core = round(Int, DATA_PACKS.core.bytes / 1e6)
    print("GREB input dataset not installed. Download whole dataset ($total MB)?
",
        "  y: all five packs ($total MB)
",
        "  n: core pack only ($core MB); optional pack prompted when a run needs them
[y/n] ")
    answer = readline(stdin)
    return lowercase(strip(answer)) in ("y", "yes")
end

"""
    download_greb_pack(pack::Symbol; dir = greb_data_dir()) -> String

Add dataset pack `pack` (`:era`, `:cmip5`, `:enso` or `:solar`, or `:all`) to the dataset directory `dir` 
and return `dir`. Downloadeds via DataDeps on first use, copies files into `dir`; skips files already present. 
An unknown pack is an `ArgumentError`.

| Pack | Adds |
|:-----|:-----|
| `:era` | ERA-Interim climatology ()`load_climatology(dir; dataset = :era)`) |
| `:cmip5` | CMIP5 RCP8.5 anomalies (`:rcp85`) |
| `:enso` | El Nino and La Nina anomalies (`:elnino`, `:lanina`) |
| `:solar` | paleo/orbital solar tables (`:paleo_231kyr`, `:eccentricity`, `:obliquity`) |
"""
function download_greb_pack(pack::Symbol; dir::Union{Nothing,AbstractString}=nothing)
    pack === :all || pack in _OPTIONAL_PACKS ||
        throw(ArgumentError("unknown dataset pack :$pack; the optional packs are " *
                            join(repr.(_OPTIONAL_PACKS), ", ") * ", or :all"))
    target = greb_data_dir(dir; packs=:core)
    if pack === :all
        foreach(p -> download_greb_pack(p; dir=target), _OPTIONAL_PACKS)
        return target
    end
    files = getfield(_OPTIONAL_PACK_FILES, pack)
    all(f -> isfile(joinpath(target, f)), files) && return target

    # Resolved lazily, like the core pack: this is what triggers the download
    source = DataDeps.resolve(_pack_dep_name(pack), @__FILE__)
    for file in files
        dest = joinpath(target, file)
        isfile(dest) && continue
        isfile(joinpath(source, file)) ||
            error("the downloaded pack :$pack lacks $file (looked in $source)")
        mkpath(dirname(dest))
        cp(joinpath(source, file), dest)
    end
    return target
end

"""
    _cached_datadep_path() -> String or nothing

Path of the already-downloaded dataset, or `nothing` if it is absent. Looks
only, never downloads: the public DataDeps calls fetch missing data, so this
uses the internal `DataDeps.try_determine_load_path`, and treats any error
from it as "not cached".
"""
function _cached_datadep_path()
    try
        path = DataDeps.try_determine_load_path(DATA_DEP_NAME, @__FILE__)
        return (path !== nothing && isdir(path)) ? String(path) : nothing
    catch
        return nothing
    end
end

"""
    register_greb_datadep()

Register the dataset packs with DataDeps. Called from `GREBClimate.__init__`;
registration itself is cheap and downloads nothing.
"""
function register_greb_datadep()
    for pack in keys(DATA_PACKS)
        register(DataDep(_pack_dep_name(pack), _pack_message(pack), _pack_url(pack), DATA_PACKS[pack].sha256;
            post_fetch_method=unpack))
    end
    return nothing
end
