# JLD2 cannot open one file from several tasks at once: runs side by side
# take this lock around every open.
const _JLD2_LOCK = ReentrantLock()

# The dataset files the loaders read, without the ".jld2" ending, keyed by the
# `ClimateFields` array each one fills. The dataset tools keep their own list
# (`tools/dataset/fields.jl`); a test checks that it holds every file named here.
const _STATIC_FILES = (z_topo="global.topography", glacier="greb.glaciers")
const _CLIMATOLOGY_FILES = (
    ncep=(Ts_clim="ncep.tsurf.1948-2007.clim", u_clim="ncep.zonal_wind.850hpa.clim",
          v_clim="ncep.meridional_wind.850hpa.clim", q_clim="ncep.atmospheric_humidity.clim",
          soil_wetness_clim="ncep.soil_moisture.clim"),
    # ERA-Interim has no soil moisture file; it uses the NCEP one
    era=(Ts_clim="erainterim.tsurf.1979-2015.clim", u_clim="erainterim.zonal_wind.850hpa.clim",
         v_clim="erainterim.meridional_wind.850hpa.clim", q_clim="erainterim.atmospheric_humidity.clim",
         soil_wetness_clim="ncep.soil_moisture.clim"),
)
const _COMMON_CLIMATOLOGY_FILES = (
    cloud_clim="isccp.cloud_cover.clim", mld_clim="woce.ocean_mixed_layer_depth.clim",
    To_clim="Tocean.clim", omega_clim="erainterim.omega.vertmean.clim",
    omega_std_clim="erainterim.omega_std.vertmean.clim", wind_speed_clim="erainterim.windspeed.850hpa.clim",
)
const _SOLAR_FILE = "solar_radiation.clim"
const _CC_ANOMALY_FILES = (
    Ts_clim_anom_cc="cmip5.tsurf.rcp85.ensmean.forcing", u_clim_anom_cc="cmip5.zonal.wind.rcp85.ensmean.forcing",
    v_clim_anom_cc="cmip5.meridional.wind.rcp85.ensmean.forcing", wind_speed_clim_anom_cc="cmip5.windspeed.rcp85.ensmean.forcing",
    omega_clim_anom_cc="cmip5.omega.rcp85.ensmean.forcing",
)
const _ENSO_EVENTS = (:elnino, :lanina)
_enso_anomaly_files(event::Symbol) = (
    Ts_clim_anom_enso="erainterim.tsurf.$event.forcing", u_clim_anom_enso="erainterim.zonal.wind.$event.forcing",
    v_clim_anom_enso="erainterim.meridional.wind.$event.forcing", wind_speed_clim_anom_enso="erainterim.windspeed.$event.forcing",
    omega_clim_anom_enso="erainterim.omega.$event.forcing",
)
# The combined file `climatology/flux_corrections.jld2`: its keys, by the array each fills
const _FLUX_CORRECTION_KEYS = (Ts_flux_correction="Tsurf_flux_correction", q_flux_correction="vapour_flux_correction",
                               To_flux_correction="Tocean_flux_correction")
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
        getfield(fields, field) .= read_field(joinpath(dir, name * ".jld2")).data
    end
end

"""
    read_field(filepath::String)

Read a `.jld2` field file written by `tools/dataset/convert_greb_to_jld2.jl`.

# Returns
- named tuple `(data, dim_names, coords, ctl)` where:
  - `data`: Array{Float32} with shape as stored
  - `dim_names`: Vector{String} of dimension names (e.g., ["lon", "lat", "time"])
  - `coords`: `Dict{Int,Vector{Float64}}` of physical coordinate values per
    dimension index, or `nothing` if the file has none
  - `ctl`: raw GrADS `.ctl` metadata text, or `nothing` if the file has none
"""
function read_field(filepath::String)
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
    load_solar_forcing(jld2_dir::String, forcing_type::Symbol, index::Int=0)

Loads an alternate solar-forcing table for paleo/orbital experiments.
`forcing_type` is `:paleo`, `:eccentricity`, or `:obliquity`; for the latter
two, `index` selects the matching row by coordinate value. Used by
[`resolve`](@ref) for a [`SolarTable`](@ref) scenario; `greb_model!` swaps it
into `fields.sw_solar` for the scenario run.
"""
function load_solar_forcing(jld2_dir::String, forcing_type::Symbol, index::Int=0)

    if forcing_type == :paleo
        filepath = joinpath(jld2_dir, "solar_scenarios", "solar_paleo.jld2")
        result = read_field(filepath)
        return result.data

    elseif forcing_type == :eccentricity
        filepath = joinpath(jld2_dir, "solar_scenarios", "solar_eccentricity.jld2")
        result = read_field(filepath)
        values = Int.(result.coords[1])
        pos = findfirst(==(index), values)
        pos === nothing && throw(ArgumentError("eccentricity index $index is not in the table; available: $(values)"))
        return result.data[pos, :, :]

    elseif forcing_type == :obliquity
        filepath = joinpath(jld2_dir, "solar_scenarios", "solar_obliquity.jld2")
        result = read_field(filepath)
        values = Int.(result.coords[1])
        pos = findfirst(==(index), values)
        pos === nothing && throw(ArgumentError("obliquity index $index is not in the table; available: $(values)"))
        return result.data[pos, :, :]

    else
        error("Unknown forcing type: $forcing_type. Use :paleo, :eccentricity, or :obliquity")
    end
end

"""
    load_co2_scenario(jld2_dir::String, scenario::Symbol) -> Dict{Int,Float32}

Loads a `year => CO2` (ppm-equivalent) lookup table for an IPCC scenario
(e.g. `:ssp585`, `:rcp85`) from the combined `scenario/ipcc_scenarios.jld2`.
The RCP6.0 table is `:rcp60`.
"""
function load_co2_scenario(jld2_dir::String, scenario::Symbol)
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
    load_co2_custom(path::String) -> Dict{Int,Float32}

Loads a `year => CO2` lookup table for the `:custom_co2` experiment from a
plain-text file, one `year CO2` pair per line. Blank lines
and lines starting with `#` are skipped.
"""
function load_co2_custom(path::String)
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
    load_flux_corrections!(jld2_dir::String, fields::ClimateFields)

Load the flux corrections from the combined `climatology/flux_corrections.jld2`
into `fields`. A missing file or a missing table in it is an `ArgumentError`.
"""
function load_flux_corrections!(jld2_dir::String, fields::ClimateFields)
    filepath = joinpath(jld2_dir, "climatology", "flux_corrections.jld2")
    isfile(filepath) || throw(ArgumentError("flux corrections file not found: $filepath"))
    @lock _JLD2_LOCK jldopen(filepath, "r") do file
        for (field, key) in pairs(_FLUX_CORRECTION_KEYS)
            haskey(file, key) || throw(ArgumentError("$key not found in $filepath"))
            getfield(fields, field) .= file[key]
            @info "Loaded $key"
        end
    end
    return nothing
end

function _load_anomaly_fields!(fields::ClimateFields, jld2_dir::String, files::NamedTuple)
    for (field, name) in pairs(files)
        filepath = joinpath(jld2_dir, "climatology", name * ".jld2")
        isfile(filepath) ||
            error("Anomaly forcing file not found: $filepath (run tools/dataset/convert_greb_to_jld2.jl)")
        getfield(fields, field) .= read_field(filepath).data
    end
end

"""
    load_boundary_anomaly!(jld2_dir::String, fields::ClimateFields, source::Symbol)

Loads the anomaly fields of a [`BoundaryAnomaly`](@ref) scenario into `fields`.
`source` is `:cmip5_rcp85` (the CMIP5 RCP8.5 ensemble mean, into
`fields.Ts_clim_anom_cc`/`u_clim_anom_cc`/`v_clim_anom_cc`/`omega_clim_anom_cc`/
`wind_speed_clim_anom_cc`), or `:elnino` or `:lanina` (the ERA-Interim composite mean,
into `fields.*_anom_enso`). Errors on a missing file rather than defaulting to
zero.

`fields` remembers the directory and the source, and [`greb_model!`](@ref)
does not read the files again for a later run on the same `fields`, directory
and source.
"""
function load_boundary_anomaly!(jld2_dir::String, fields::ClimateFields, source::Symbol)
    if source === :cmip5_rcp85
        fields.anom_cc_source = ""
        _load_anomaly_fields!(fields, jld2_dir, _CC_ANOMALY_FILES)
        fields.anom_cc_source = jld2_dir
    elseif source in _ENSO_EVENTS
        fields.anom_enso_source = ("", :none)
        _load_anomaly_fields!(fields, jld2_dir, _enso_anomaly_files(source))
        fields.anom_enso_source = (jld2_dir, source)
    else
        throw(ArgumentError("source must be :cmip5_rcp85, :elnino or :lanina, got :$source"))
    end
    return nothing
end

"""
    load_climatology(jld2_dir::String; dataset::Symbol=:ncep, corrections::Bool=true)

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
function load_climatology(jld2_dir::String; dataset::Symbol=:ncep, corrections::Bool=true)
    if !isdir(jld2_dir)
        error("JLD2 directory not found: $jld2_dir")
    end

    fields = ClimateFields()

    haskey(_CLIMATOLOGY_FILES, dataset) ||
        throw(ArgumentError("unknown dataset :$dataset; use one of $(join(repr.(keys(_CLIMATOLOGY_FILES)), ", "))"))

    @info "Loading static fields"
    _load_fields!(fields, joinpath(jld2_dir, "static"), _STATIC_FILES)

    @info "Loading the $dataset climatology"
    climatology_dir = joinpath(jld2_dir, "climatology")
    _load_fields!(fields, climatology_dir, _CLIMATOLOGY_FILES[dataset])

    @info "Loading the common climatology fields"
    _load_fields!(fields, climatology_dir, _COMMON_CLIMATOLOGY_FILES)

    # Solar radiation (special: lat × time)
    @info "Loading the solar radiation table"
    solar_path = joinpath(jld2_dir, "solar", _SOLAR_FILE * ".jld2")
    if isfile(solar_path)
        solar_result = read_field(solar_path)
        size(solar_result.data) == (ydim, nstep_yr) ||
            error("$solar_path holds a $(size(solar_result.data)) table, expected ($ydim, $nstep_yr)")
        fields.sw_solar .= solar_result.data
    else
        error("Solar radiation file not found: $solar_path")
    end

    if corrections
        @info "Loading the flux corrections"
        load_flux_corrections!(jld2_dir, fields)
    end

    split_winds!(fields)

    fields.loaded = true
    @info "Dataset loaded from $jld2_dir"
    return fields
end
