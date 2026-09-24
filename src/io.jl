"""
    read_jld2(filepath::String)

Read a `.jld2` field file written by `tools/convert_greb_to_jld2.jl`.

# Returns
- named tuple `(data, dim_names, coords, ctl)` where:
  - `data`: Array{Float32} with shape as stored
  - `dim_names`: Vector{String} of dimension names (e.g., ["lon", "lat", "time"])
  - `coords`: `Dict{Int,Vector{Float64}}` of physical coordinate values per
    dimension index, or `nothing` if the file has none
  - `ctl`: raw GrADS `.ctl` metadata text, or `nothing` if the file has none
"""
function read_jld2(filepath::String)
    jldopen(filepath, "r") do file
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
`greb_model!` to temporarily swap `fields.sw_solar` for these experiments.
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
        @assert pos !== nothing "Eccentricity index $index not found in $(values)"
        return result.data[pos, :, :]

    elseif forcing_type == :obliquity
        filepath = joinpath(jld2_dir, "solar_scenarios", "solar_obliquity.jld2")
        result = read_jld2(filepath)
        values = Int.(result.coords[1])
        pos = findfirst(==(index), values)
        @assert pos !== nothing "Obliquity index $index not found in $(values)"
        return result.data[pos, :, :]

    else
        error("Unknown forcing type: $forcing_type. Use :paleo, :eccentricity, or :obliquity")
    end
end

"""
    load_co2_scenario_jld2(jld2_dir::String, scenario::Symbol) -> Dict{Int,Float32}

Loads a `year => CO2` (ppm-equivalent) lookup table for an IPCC scenario
(e.g. `:ssp585`, `:rcp85`) from the combined `scenario/ipcc_scenarios.jld2`.
"""
function load_co2_scenario_jld2(jld2_dir::String, scenario::Symbol)
    filepath = joinpath(jld2_dir, "scenario", "ipcc_scenarios.jld2")
    isfile(filepath) ||
        error("Scenario file not found: $filepath (run tools/convert_greb_to_jld2.jl)")
    scenarios = jldopen(filepath, "r") do file
        file["scenarios"]
    end
    key = string(scenario)
    haskey(scenarios, key) ||
        error("No CO2 scenario table for \"$key\" in $filepath. Available: $(sort(collect(keys(scenarios))))")
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
    _each_year_value(path, "custom CO2 scenario", "year CO2") do yr, v
        table[yr] = parse(Float32, v)
    end
    return table
end

# Calls `f(year, value)` for each `year value` line of a text table, skipping
# blank lines and `#` comments; `value` is the unparsed second column.
function _each_year_value(f, path, label, columns)
    open(path) do io
        for line in eachline(io)
            stripped = strip(line)
            (isempty(stripped) || startswith(stripped, "#")) && continue
            cols = split(stripped)
            length(cols) >= 2 ||
                error("Malformed line in $label file $path: \"$line\" (expected \"$columns\")")
            f(parse(Int, cols[1]), cols[2])
        end
    end
end

"""
    load_solar_series(path; reference=nothing) -> Dict{Int,Float32}

Load a `year TSI` table from a plain-text file (blank lines and `#` skipped)
for `cfg.solar_scenario`. Values are divided by `reference`, defaulting to the
table mean, so the multiplier averages 1 and the control climate is unchanged.
Pass `reference` (W/m2) to normalise against a solar constant instead.
Composes with any experiment; does not replace the built-in solar experiments.
"""
function load_solar_series(path::AbstractString; reference::Union{Real,Nothing}=nothing)
    isfile(path) || error("Solar series file not found: $path")
    raw = Dict{Int,Float64}()
    _each_year_value(path, "solar series", "year TSI") do yr, v
        haskey(raw, yr) && error("Solar series file $path lists year $yr more than once")
        raw[yr] = parse(Float64, v)
    end
    isempty(raw) && error("Solar series file $path contains no data lines")
    # A missing year would otherwise surface as an error in the middle of a run.
    yrs = sort!(collect(keys(raw)))
    for (a, b) in zip(yrs, yrs[2:end])
        b == a + 1 || error("Solar series $path has a gap: $a is followed by $b")
    end
    ref = reference === nothing ? sum(values(raw)) / length(raw) : Float64(reference)
    ref > 0 || error("Solar series reference must be positive, got $ref")
    return Dict{Int,Float32}(yr => Float32(v / ref) for (yr, v) in raw)
end

# Forcing-series standard: the one JLD2 layout for time-varying inputs (spec in
# the vault, 07-standards/forcing-series-format). Times are decimal years on
# the model's own calendar, so the model never converts dates.
const FORCING_SERIES_FORMAT = 2
const FORCING_CALENDAR = "greb_365"

# Reads and validates a forcing-series file of `kind` with `units`. Returns
# `(time, lat, values, source)`; `lat` is `nothing` for a global series, in
# which case `values` is a vector over time instead of a latitude-by-time matrix.
function _read_forcing_series(path, kind, units)
    isfile(path) || error("Forcing series file not found: $path")
    fix = "regenerate it with the converter in tools/"
    time, lat, values, source = jldopen(path) do file
        for key in ("format_version", "kind", "calendar", "time", "values", "units")
            haskey(file, key) || error("Forcing series file $path is missing `$key`; $fix")
        end
        v = file["format_version"]
        v == FORCING_SERIES_FORMAT ||
            error("Forcing series file $path has format_version $v; this version reads " *
                  "$FORCING_SERIES_FORMAT; $fix")
        file["kind"] == kind || error("Forcing series file $path holds `$(file["kind"])`, not `$kind`")
        file["calendar"] == FORCING_CALENDAR ||
            error("Forcing series file $path uses calendar `$(file["calendar"])`, not `$FORCING_CALENDAR`")
        file["units"] == units ||
            error("Forcing series file $path has units `$(file["units"])`, expected `$units`")
        (Float64.(file["time"]), haskey(file, "lat") ? Float64.(file["lat"]) : nothing,
         Float64.(file["values"]), haskey(file, "source") ? String(file["source"]) : "")
    end
    isempty(time) && error("Forcing series file $path contains no times")
    _strictly_ascending(time) || error("Forcing series file $path: `time` must be strictly ascending")
    expected = lat === nothing ? (length(time),) : (length(lat), length(time))
    size(values) == expected ||
        error("Forcing series file $path: `values` has size $(size(values)), expected $expected")
    lat === nothing || _strictly_ascending(lat) ||
        error("Forcing series file $path: `lat` must be strictly ascending")
    all(isfinite, values) || error("Forcing series file $path: `values` contains non-finite entries")
    return time, lat, values, source
end

"""
    load_aerosol_series(path::AbstractString) -> AerosolSeries

Loads a stratospheric optical-depth record for `AerosolScenario(series=...)`
from a forcing-series file of kind `aerosol_aod` (units `1`, per latitude),
as written by `tools/convert_aerosol_to_jld2.jl`. Values are interpolated
linearly onto the model latitudes (end values held) and, during the run, in
time. The record is zero outside its first and last time and is never
altered, so loading warns when it ends above `SERIES_END_WARN_AOD`: the
dimming would stop abruptly.
"""
function load_aerosol_series(path::AbstractString)
    years, lat, aod, source = _read_forcing_series(path, "aerosol_aod", "1")
    lat === nothing && error("Aerosol series file $path has no `lat`; aerosol records are zonal")
    out = Matrix{Float32}(undef, ydim, length(years))
    for k in eachindex(years), j in 1:ydim
        out[j, k] = Float32(_interp_lat(lat, @view(aod[:, k]), Float64(lat_grid[j])))
    end
    last_aod = maximum(@view aod[:, end])
    last_aod > SERIES_END_WARN_AOD &&
        @warn "Aerosol series $path ends at year $(years[end]) with optical depth " *
              "$last_aod; the record is zero after its last year, so the dimming stops abruptly"
    return AerosolSeries(years, out, source)
end

_strictly_ascending(v) = all(i -> v[i] < v[i + 1], firstindex(v):lastindex(v) - 1)

# Linear interpolation in `xs` (ascending), holding the end values beyond it.
function _interp_lat(xs, ys, x)
    x <= first(xs) && return first(ys)
    x >= last(xs) && return last(ys)
    k = searchsortedlast(xs, x)
    f = (x - xs[k]) / (xs[k + 1] - xs[k])
    return (1 - f) * ys[k] + f * ys[k + 1]
end

"""
Load flux corrections from the combined `climatology/flux_corrections.jld2`
into `fields` (zeros per-field if the file or an individual key is missing).
"""
function load_flux_corrections_jld2!(jld2_dir::String, fields::ClimateFields)
    correction_keys = Dict(
        "Tsurf_flux_correction" => fields.TF_correct,
        "vapour_flux_correction" => fields.qF_correct,
        "Tocean_flux_correction" => fields.ToF_correct
    )

    filepath = joinpath(jld2_dir, "climatology", "flux_corrections.jld2")
    if isfile(filepath)
        jldopen(filepath, "r") do file
            for (key, array) in correction_keys
                if haskey(file, key)
                    array .= file[key]
                    println("✅ Loaded $key")
                else
                    fill!(array, 0.0)
                    @warn "$key not found in $filepath, using zeros"
                end
            end
        end
    else
        for array in values(correction_keys)
            fill!(array, 0.0)
        end
        @warn "$filepath not found, using zeros for all flux corrections"
    end
end

function _load_anomaly_field!(climatology_dir::String, filename::String, target::AbstractArray)
    filepath = joinpath(climatology_dir, filename)
    isfile(filepath) ||
        error("Anomaly forcing file not found: $filepath (run tools/convert_greb_to_jld2.jl)")
    target .= read_jld2(filepath).data
end

"""
    load_cc_anomaly_jld2!(jld2_dir::String, fields::ClimateFields, cfg::PhysicsConfig)

Loads the CMIP5 RCP8.5 ensemble-mean climate-change anomaly fields into
`fields.Tclim_anom_cc`/`uclim_anom_cc`/`vclim_anom_cc`/`omegaclim_anom_cc`/
`wsclim_anom_cc`, gated per-field by `cfg.log_tsurf_ext`/`log_hwind_ext`/
`log_omega_ext`. Errors on a missing file rather than defaulting to zero, since this data
*is* the `:rcp85` experiment's forcing.
"""
function load_cc_anomaly_jld2!(jld2_dir::String, fields::ClimateFields, cfg::PhysicsConfig)
    dir = joinpath(jld2_dir, "climatology")
    cfg.log_tsurf_ext && _load_anomaly_field!(dir, "cmip5.tsurf.rcp85.ensmean.forcing.jld2", fields.Tclim_anom_cc)
    if cfg.log_hwind_ext
        _load_anomaly_field!(dir, "cmip5.zonal.wind.rcp85.ensmean.forcing.jld2", fields.uclim_anom_cc)
        _load_anomaly_field!(dir, "cmip5.meridional.wind.rcp85.ensmean.forcing.jld2", fields.vclim_anom_cc)
        _load_anomaly_field!(dir, "cmip5.windspeed.rcp85.ensmean.forcing.jld2", fields.wsclim_anom_cc)
    end
    cfg.log_omega_ext && _load_anomaly_field!(dir, "cmip5.omega.rcp85.ensmean.forcing.jld2", fields.omegaclim_anom_cc)
end

"""
    load_enso_anomaly_jld2!(jld2_dir::String, fields::ClimateFields, cfg::PhysicsConfig, which::Symbol)

Loads the ERA-Interim composite-mean El Niño (`which=:elnino`) or La Niña
(`:lanina`) anomaly fields into `fields.*_anom_enso`, gated the same way as
[`load_cc_anomaly_jld2!`](@ref).
"""
function load_enso_anomaly_jld2!(jld2_dir::String, fields::ClimateFields, cfg::PhysicsConfig, which::Symbol)
    suffix = which == :elnino ? "elnino" :
             which == :lanina ? "lanina" :
             error("which must be :elnino or :lanina, got $which")
    dir = joinpath(jld2_dir, "climatology")
    cfg.log_tsurf_ext && _load_anomaly_field!(dir, "erainterim.tsurf.$suffix.forcing.jld2", fields.Tclim_anom_enso)
    if cfg.log_hwind_ext
        _load_anomaly_field!(dir, "erainterim.zonal.wind.$suffix.forcing.jld2", fields.uclim_anom_enso)
        _load_anomaly_field!(dir, "erainterim.meridional.wind.$suffix.forcing.jld2", fields.vclim_anom_enso)
        _load_anomaly_field!(dir, "erainterim.windspeed.$suffix.forcing.jld2", fields.wsclim_anom_enso)
    end
    cfg.log_omega_ext && _load_anomaly_field!(dir, "erainterim.omega.$suffix.forcing.jld2", fields.omegaclim_anom_enso)
end

"""
    load_greb_jld2!(jld2_dir::String; dataset::Symbol=:ncep)

Load all GREB input data from JLD2 formatted files, returning a fresh
[`ClimateFields`](@ref). `dataset` (`:ncep`/`:era`) selects which
climatology *files* to read; this is independent of `PhysicsConfig.log_clim`,
which only selects hydrology regression *coefficients* in
[`set_hydrology_parameters!`](@ref).
"""
function load_greb_jld2!(jld2_dir::String; dataset::Symbol=:ncep)
    if !isdir(jld2_dir)
        error("JLD2 directory not found: $jld2_dir")
    end

    fields = ClimateFields()

    # Static 2D files
    println("📂 Loading static fields...")
    topo_result = read_jld2(joinpath(jld2_dir, "static", "global.topography.jld2"))
    fields.z_topo .= topo_result.data

    glacier_result = read_jld2(joinpath(jld2_dir, "static", "greb.glaciers.jld2"))
    fields.glacier .= glacier_result.data

    # Dataset-specific file mapping
    file_map = Dict(
        :ncep => Dict(
            "Tclim" => "ncep.tsurf.1948-2007.clim.jld2",
            "uclim" => "ncep.zonal_wind.850hpa.clim.jld2",
            "vclim" => "ncep.meridional_wind.850hpa.clim.jld2",
            "qclim" => "ncep.atmospheric_humidity.clim.jld2",
            "swetclim" => "ncep.soil_moisture.clim.jld2"
        ),
        :era => Dict(
            "Tclim" => "erainterim.tsurf.1979-2015.clim.jld2",
            "uclim" => "erainterim.zonal_wind.850hpa.clim.jld2",
            "vclim" => "erainterim.meridional_wind.850hpa.clim.jld2",
            "qclim" => "erainterim.atmospheric_humidity.clim.jld2",
            "swetclim" => "ncep.soil_moisture.clim.jld2"
        )
    )

    # Use mixed dataset as fallback
    files = get(file_map, dataset, file_map[:ncep])

    println("📂 Loading 3D climatology ($dataset dataset)...")
    climatology_dir = joinpath(jld2_dir, "climatology")

    # Load each variable individually with unique result names
    tsurf_result = read_jld2(joinpath(climatology_dir, files["Tclim"]))
    fields.Tclim .= tsurf_result.data

    uwind_result = read_jld2(joinpath(climatology_dir, files["uclim"]))
    fields.uclim .= uwind_result.data

    vwind_result = read_jld2(joinpath(climatology_dir, files["vclim"]))
    fields.vclim .= vwind_result.data

    humid_result = read_jld2(joinpath(climatology_dir, files["qclim"]))
    fields.qclim .= humid_result.data

    swet_result = read_jld2(joinpath(climatology_dir, files["swetclim"]))
    fields.swetclim .= swet_result.data

    # Common climatology files
    println("📂 Loading common climatology fields...")

    cld_result = read_jld2(joinpath(climatology_dir, "isccp.cloud_cover.clim.jld2"))
    fields.cldclim .= cld_result.data

    mld_result = read_jld2(joinpath(climatology_dir, "woce.ocean_mixed_layer_depth.clim.jld2"))
    fields.mldclim .= mld_result.data

    tocean_result = read_jld2(joinpath(climatology_dir, "Tocean.clim.jld2"))
    fields.Toclim .= tocean_result.data

    omega_result = read_jld2(joinpath(climatology_dir, "erainterim.omega.vertmean.clim.jld2"))
    fields.omegaclim .= omega_result.data

    omegastd_result = read_jld2(joinpath(climatology_dir, "erainterim.omega_std.vertmean.clim.jld2"))
    fields.omegastdclim .= omegastd_result.data

    ws_result = read_jld2(joinpath(climatology_dir, "erainterim.windspeed.850hpa.clim.jld2"))
    fields.wsclim .= ws_result.data

    # Solar radiation (special: lat × time)
    println("📂 Loading solar radiation...")
    solar_path = joinpath(jld2_dir, "solar", "solar_radiation.clim.jld2")
    if isfile(solar_path)
        solar_result = read_jld2(solar_path)
        @assert size(solar_result.data) == (ydim, nstep_yr) "Wrong solar dimensions"
        fields.sw_solar .= solar_result.data
    else
        error("Solar radiation file not found: $solar_path")
    end

    # Optional: Load flux corrections
    println("📂 Loading flux corrections...")
    load_flux_corrections_jld2!(jld2_dir, fields)

    # Update wind sign splits
    @. fields.uclim_m = ifelse(fields.uclim >= 0.0, fields.uclim, 0.0)
    @. fields.uclim_p = ifelse(fields.uclim < 0.0, fields.uclim, 0.0)
    @. fields.vclim_m = ifelse(fields.vclim >= 0.0, fields.vclim, 0.0)
    @. fields.vclim_p = ifelse(fields.vclim < 0.0, fields.vclim, 0.0)

    fields.loaded = true
    println("✅ All GREB data loaded successfully from JLD2")
    return fields
end
