# JLD2 cannot open one file from several tasks at once: runs side by side
# take this lock around every open.
const _JLD2_LOCK = ReentrantLock()

# The dataset files the loaders read, without the ".jld2" ending, keyed by the
# `ClimateFields` array each one fills. The dataset tools keep their own list
# (`tools/dataset/fields.jl`); a test checks that it holds every file named here.
const _STATIC_FILES = (z_topo="global.topography", glacier="greb.glaciers")
const _CLIMATOLOGY_FILES = (
    ncep=(Tclim="ncep.tsurf.1948-2007.clim", uclim="ncep.zonal_wind.850hpa.clim",
          vclim="ncep.meridional_wind.850hpa.clim", qclim="ncep.atmospheric_humidity.clim",
          swetclim="ncep.soil_moisture.clim"),
    # ERA-Interim has no soil moisture file; it uses the NCEP one
    era=(Tclim="erainterim.tsurf.1979-2015.clim", uclim="erainterim.zonal_wind.850hpa.clim",
         vclim="erainterim.meridional_wind.850hpa.clim", qclim="erainterim.atmospheric_humidity.clim",
         swetclim="ncep.soil_moisture.clim"),
)
const _COMMON_CLIMATOLOGY_FILES = (
    cldclim="isccp.cloud_cover.clim", mldclim="woce.ocean_mixed_layer_depth.clim",
    Toclim="Tocean.clim", omegaclim="erainterim.omega.vertmean.clim",
    omegastdclim="erainterim.omega_std.vertmean.clim", wsclim="erainterim.windspeed.850hpa.clim",
)
const _SOLAR_FILE = "solar_radiation.clim"
const _CC_ANOMALY_FILES = (
    Tclim_anom_cc="cmip5.tsurf.rcp85.ensmean.forcing", uclim_anom_cc="cmip5.zonal.wind.rcp85.ensmean.forcing",
    vclim_anom_cc="cmip5.meridional.wind.rcp85.ensmean.forcing", wsclim_anom_cc="cmip5.windspeed.rcp85.ensmean.forcing",
    omegaclim_anom_cc="cmip5.omega.rcp85.ensmean.forcing",
)
const _ENSO_EVENTS = (:elnino, :lanina)
_enso_anomaly_files(event::Symbol) = (
    Tclim_anom_enso="erainterim.tsurf.$event.forcing", uclim_anom_enso="erainterim.zonal.wind.$event.forcing",
    vclim_anom_enso="erainterim.meridional.wind.$event.forcing", wsclim_anom_enso="erainterim.windspeed.$event.forcing",
    omegaclim_anom_enso="erainterim.omega.$event.forcing",
)
# The combined file `climatology/flux_corrections.jld2`: its keys, by the array each fills
const _FLUX_CORRECTION_KEYS = (TF_correct="Tsurf_flux_correction", qF_correct="vapour_flux_correction",
                               ToF_correct="Tocean_flux_correction")
# Files that hold several tables each and are read by their own loaders
const _COMBINED_FILES = ("flux_corrections", "ipcc_scenarios", "solar_paleo", "solar_eccentricity", "solar_obliquity")

"""
    dataset_field_files() -> Set{String}

Names (without `.jld2`) of the single-field files the loaders read: the two
static fields, the solar table, the climatologies of both datasets and the
anomaly files. The dataset tools must convert and package at least these.
"""
function dataset_field_files()
    names = Set{String}(values(_STATIC_FILES))
    push!(names, _SOLAR_FILE)
    for files in (values(_CLIMATOLOGY_FILES)..., _COMMON_CLIMATOLOGY_FILES, _CC_ANOMALY_FILES,
                  (_enso_anomaly_files(e) for e in _ENSO_EVENTS)...)
        union!(names, values(files))
    end
    return names
end

# Fill the arrays of `fields` named in `files` from the files in `dir`
function _load_fields!(fields::ClimateFields, dir::String, files::NamedTuple)
    for (field, name) in pairs(files)
        getfield(fields, field) .= read_jld2(joinpath(dir, name * ".jld2")).data
    end
end

"""
    read_jld2(filepath::String)

Read a `.jld2` field file written by `tools/dataset/convert_greb_to_jld2.jl`.

# Returns
- named tuple `(data, dim_names, coords, ctl)` where:
  - `data`: Array{Float32} with shape as stored
  - `dim_names`: Vector{String} of dimension names (e.g., ["lon", "lat", "time"])
  - `coords`: `Dict{Int,Vector{Float64}}` of physical coordinate values per
    dimension index, or `nothing` if the file has none
  - `ctl`: raw GrADS `.ctl` metadata text, or `nothing` if the file has none
"""
function read_jld2(filepath::String)
    @lock _JLD2_LOCK jldopen(filepath, "r") do file
        return (
            data=file["data"],
            dim_names=file["dim_names"],
            coords=haskey(file, "coords") ? file["coords"] : nothing,
            ctl=haskey(file, "ctl") ? file["ctl"] : nothing,
        )
    end
end

"""
    load_solar_forcing_jld2(jld2_dir::String, forcing_type::Symbol, index::Int=0)

Loads an alternate solar-forcing table for paleo/orbital experiments.
`forcing_type` is `:paleo`, `:eccentricity`, or `:obliquity`; for the latter
two, `index` selects the matching row by coordinate value. Used by
[`resolve`](@ref) for a [`SolarTable`](@ref) scenario; `greb_model!` swaps it
into `fields.sw_solar` for the scenario run.
"""
function load_solar_forcing_jld2(jld2_dir::String, forcing_type::Symbol, index::Int=0)

    if forcing_type == :paleo
        filepath = joinpath(jld2_dir, "solar_scenarios", "solar_paleo.jld2")
        result = read_jld2(filepath)
        return result.data

    elseif forcing_type == :eccentricity
        filepath = joinpath(jld2_dir, "solar_scenarios", "solar_eccentricity.jld2")
        result = read_jld2(filepath)
        values = Int.(result.coords[1])
        pos = findfirst(==(index), values)
        pos === nothing && throw(ArgumentError("eccentricity index $index is not in the table; available: $(values)"))
        return result.data[pos, :, :]

    elseif forcing_type == :obliquity
        filepath = joinpath(jld2_dir, "solar_scenarios", "solar_obliquity.jld2")
        result = read_jld2(filepath)
        values = Int.(result.coords[1])
        pos = findfirst(==(index), values)
        pos === nothing && throw(ArgumentError("obliquity index $index is not in the table; available: $(values)"))
        return result.data[pos, :, :]

    else
        error("Unknown forcing type: $forcing_type. Use :paleo, :eccentricity, or :obliquity")
    end
end

"""
    load_co2_scenario_jld2(jld2_dir::String, scenario::Symbol) -> Dict{Int,Float32}

Loads a `year => CO2` (ppm-equivalent) lookup table for an IPCC scenario
(e.g. `:ssp585`, `:rcp85`) from the combined `scenario/ipcc_scenarios.jld2`.
The RCP6.0 table is `:rcp60`.
"""
function load_co2_scenario_jld2(jld2_dir::String, scenario::Symbol)
    filepath = joinpath(jld2_dir, "scenario", "ipcc_scenarios.jld2")
    isfile(filepath) ||
        error("Scenario file not found: $filepath (run tools/dataset/convert_greb_to_jld2.jl)")
    scenarios = @lock _JLD2_LOCK jldopen(filepath, "r") do file
        file["scenarios"]
    end
    # The dataset stores the RCP6.0 table as "rcp6"; its name here is :rcp60
    scenario === :rcp6 && throw(ArgumentError("the RCP6.0 table is :rcp60, not :rcp6"))
    key = scenario === :rcp60 ? "rcp6" : string(scenario)
    haskey(scenarios, key) ||
        error("No CO2 scenario table for :$scenario in $filepath. Available: $(sort(collect(keys(scenarios))))")
    return Dict{Int,Float32}(yr => Float32(co2) for (yr, co2) in scenarios[key])
end

"""
    load_custom_co2_scenario(path::String) -> Dict{Int,Float32}

Loads a `year => CO2` lookup table for the `:custom_co2` experiment from a
plain-text file, one `year CO2` pair per line. Blank lines
and lines starting with `#` are skipped.
"""
function load_custom_co2_scenario(path::String)
    isfile(path) || error("Custom CO2 scenario file not found: $path")
    table = Dict{Int,Float32}()
    open(path) do io
        for line in eachline(io)
            stripped = strip(line)
            (isempty(stripped) || startswith(stripped, "#")) && continue
            cols = split(stripped)
            length(cols) >= 2 ||
                error("Malformed line in custom CO2 scenario file $path: \"$line\" (expected \"year CO2\")")
            table[parse(Int, cols[1])] = parse(Float32, cols[2])
        end
    end
    return table
end

"""
    load_flux_corrections_jld2!(jld2_dir::String, fields::ClimateFields)

Load the flux corrections from the combined `climatology/flux_corrections.jld2`
into `fields`. A missing file or a missing table in it is an `ArgumentError`.
"""
function load_flux_corrections_jld2!(jld2_dir::String, fields::ClimateFields)
    filepath = joinpath(jld2_dir, "climatology", "flux_corrections.jld2")
    isfile(filepath) || throw(ArgumentError("flux corrections file not found: $filepath"))
    @lock _JLD2_LOCK jldopen(filepath, "r") do file
        for (field, key) in pairs(_FLUX_CORRECTION_KEYS)
            haskey(file, key) || throw(ArgumentError("$key not found in $filepath"))
            getfield(fields, field) .= file[key]
            println("✅ Loaded $key")
        end
    end
    return nothing
end

function _load_anomaly_fields!(fields::ClimateFields, jld2_dir::String, files::NamedTuple)
    for (field, name) in pairs(files)
        filepath = joinpath(jld2_dir, "climatology", name * ".jld2")
        isfile(filepath) ||
            error("Anomaly forcing file not found: $filepath (run tools/dataset/convert_greb_to_jld2.jl)")
        getfield(fields, field) .= read_jld2(filepath).data
    end
end

"""
    load_cc_anomaly_jld2!(jld2_dir::String, fields::ClimateFields)

Loads the CMIP5 RCP8.5 ensemble-mean climate-change anomaly fields into
`fields.Tclim_anom_cc`/`uclim_anom_cc`/`vclim_anom_cc`/`omegaclim_anom_cc`/
`wsclim_anom_cc`: the forcing of a `BoundaryAnomaly(:cmip5_rcp85)` scenario.
Errors on a missing file rather than defaulting to zero.

`fields` remembers the directory, and [`greb_model!`](@ref) does not read the
files again for a later run on the same `fields` and directory.
"""
function load_cc_anomaly_jld2!(jld2_dir::String, fields::ClimateFields)
    fields.anom_cc_source = ""
    _load_anomaly_fields!(fields, jld2_dir, _CC_ANOMALY_FILES)
    fields.anom_cc_source = jld2_dir
    return nothing
end

"""
    load_enso_anomaly_jld2!(jld2_dir::String, fields::ClimateFields, which::Symbol)

Loads the ERA-Interim composite-mean El Niño (`which=:elnino`) or La Niña
(`:lanina`) anomaly fields into `fields.*_anom_enso`, as
[`load_cc_anomaly_jld2!`](@ref) does for RCP8.5; `fields` remembers the
directory and the event.
"""
function load_enso_anomaly_jld2!(jld2_dir::String, fields::ClimateFields, which::Symbol)
    which in _ENSO_EVENTS || error("which must be :elnino or :lanina, got $which")
    fields.anom_enso_source = ("", :none)
    _load_anomaly_fields!(fields, jld2_dir, _enso_anomaly_files(which))
    fields.anom_enso_source = (jld2_dir, which)
    return nothing
end

"""
    load_greb_jld2!(jld2_dir::String; dataset::Symbol=:ncep, corrections::Bool=true)

Load all GREB input data from JLD2 formatted files, returning a fresh
[`ClimateFields`](@ref). `dataset` (`:ncep`/`:era`) selects which
climatology *files* to read, and any other value is an `ArgumentError`; this
is independent of `Hydrology.rain_fit`, which only selects the
rain-regression *coefficients* (see [`Hydrology`](@ref)).

Every file must be there: a missing one is an error. `corrections = false`
leaves out the stored flux corrections (`climatology/flux_corrections.jld2`),
for a dataset that has none: the three correction arrays stay zero, and a run
computes its own with [`SpinUp`](@ref) or runs without, with
[`NoCorrections`](@ref).
"""
function load_greb_jld2!(jld2_dir::String; dataset::Symbol=:ncep, corrections::Bool=true)
    if !isdir(jld2_dir)
        error("JLD2 directory not found: $jld2_dir")
    end

    fields = ClimateFields()

    haskey(_CLIMATOLOGY_FILES, dataset) ||
        throw(ArgumentError("unknown dataset :$dataset; use one of $(join(repr.(keys(_CLIMATOLOGY_FILES)), ", "))"))

    println("📂 Loading static fields...")
    _load_fields!(fields, joinpath(jld2_dir, "static"), _STATIC_FILES)

    println("📂 Loading 3D climatology ($dataset dataset)...")
    climatology_dir = joinpath(jld2_dir, "climatology")
    _load_fields!(fields, climatology_dir, _CLIMATOLOGY_FILES[dataset])

    println("📂 Loading common climatology fields...")
    _load_fields!(fields, climatology_dir, _COMMON_CLIMATOLOGY_FILES)

    # Solar radiation (special: lat × time)
    println("📂 Loading solar radiation...")
    solar_path = joinpath(jld2_dir, "solar", _SOLAR_FILE * ".jld2")
    if isfile(solar_path)
        solar_result = read_jld2(solar_path)
        size(solar_result.data) == (ydim, nstep_yr) ||
            error("$solar_path holds a $(size(solar_result.data)) table, expected ($ydim, $nstep_yr)")
        fields.sw_solar .= solar_result.data
    else
        error("Solar radiation file not found: $solar_path")
    end

    if corrections
        println("📂 Loading flux corrections...")
        load_flux_corrections_jld2!(jld2_dir, fields)
    end

    split_winds!(fields)

    fields.loaded = true
    println("✅ All GREB data loaded successfully from JLD2")
    return fields
end
