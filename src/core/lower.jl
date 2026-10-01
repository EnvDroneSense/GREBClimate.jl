# Translation of a v2 configuration into the legacy PhysicsConfig, which the
# model still runs on. Temporary: it goes when the kernels read the v2 parts.

const _RAIN_CODE = Dict(:original => -1, :fitted => 0, :rh => 1, :omega => 2, :rh_omega => 3)
const _EVAPORATION_CODE = Dict(:original => -1, :skin => 0, :original_gust => 1, :skin_gust => 2)

function _lower_hydrology!(cfg::PhysicsConfig, h::Hydrology)
    cfg.log_rain = _RAIN_CODE[h.rain]
    cfg.log_eva = _EVAPORATION_CODE[h.evaporation]
    cfg.log_clim = h.rain_fit === :ncep ? 1 : 0
    return cfg
end

function _lower_corrections!(cfg::PhysicsConfig, c::Corrections)
    cfg.flux_corrections = c isa SpinUp ? :spinup : c isa Stored ? :stored : :none
    return cfg
end

function _lower_switches!(cfg::PhysicsConfig, p::Processes, c::Corrections)
    cfg.log_atmos_dmc = p.atmosphere
    cfg.log_clouds_dmc = p.clouds !== :none
    cfg.log_clouds_drsp = p.clouds !== :uniform
    cfg.log_humid_drsp = p.humidity !== :uniform
    cfg.log_hydro_dmc = p.hydrology !== :none
    cfg.log_hydro_drsp = p.hydrology !== :no_evap_rain
    cfg.log_ocean_dmc = p.ocean !== :none
    cfg.log_ocean_drsp = p.ocean !== :mixed_layer
    cfg.log_topo_drsp = p.topography === :observed
    cfg.log_co2_dmc = p.co2
    cfg.log_ice = p.ice_albedo
    cfg.log_crcl_dmc = cfg.log_crcl_drsp = p.transport
    cfg.log_hdif = p.heat_diffusion
    cfg.log_hadv = p.heat_advection
    cfg.log_vdif = p.vapour_diffusion
    cfg.log_vadv = p.vapour_advection
    cfg.log_conv = p.moisture_convergence
    return _lower_corrections!(cfg, c)
end

# Preset scenario => legacy experiment symbol. The deconstruction presets run
# the forcing of :full_model and :co2_double.
const _LEGACY_NAME = Dict(:rcp85_boundary => :rcp85, :co2_abrupt_reverse => :co2_step, :a1b => :a1b_scenario)
const _SCENARIO_TO_LEGACY = Dict(s => get(_LEGACY_NAME, p, p)
                                 for (p, s) in _PRESET_SCENARIOS if !(p in (:decon_mean_climate, :decon_2xco2)))

# The scenario with its free parameters reset, to find its preset
function _template(s::Scenario)
    s.solar isa SolarTable && (s = _with(s; solar=SolarTable(s.solar.kind)))
    s.solar isa EarthSunDistance && (s = _with(s; solar=EarthSunDistance(0)))
    s.co2 isa CO2File && (s = _with(s; co2=CO2File("")))
    return s
end

function _lower_scenario!(cfg::PhysicsConfig, s::Scenario)
    experiment = get(_SCENARIO_TO_LEGACY, _template(s), nothing)
    experiment === nothing && throw(ArgumentError(
        "this scenario matches no preset; free combinations of its parts are not available yet"))
    cfg.experiment = experiment
    for (field, value) in pairs(_EXPERIMENT_OVERRIDES[experiment])
        setfield!(cfg, field, value)
    end
    s.solar isa SolarTable && (cfg.orbital_index = s.solar.index)
    s.solar isa EarthSunDistance && (cfg.earth_sun_distance_pct = s.solar.pct)
    s.co2 isa CO2File && (cfg.custom_co2_path = s.co2.path)
    return cfg
end

# The legacy PhysicsConfig a Config runs as
function _lower(c::Config)
    isempty(c.modules) || throw(ArgumentError("add-on modules are not available yet"))
    cfg = PhysicsConfig()
    _lower_scenario!(cfg, c.scenario)
    _lower_switches!(cfg, c.processes, c.corrections)
    return _lower_hydrology!(cfg, c.hydrology)
end

"""
    greb_model!(run::RunSpec, config::Config; jld2_dir="", fields=ClimateFields(),
                allow_uninitialized=false)

Run `config` for `run.ctrl` control and `run.scnr` scenario years. The
spin-up length comes from `config.corrections` ([`SpinUp`](@ref)`(years)`);
`run.flux` is not used. Returns `(ctrl, scnr)` as the `PhysicsConfig` method
does.
"""
function greb_model!(run::RunSpec, config::Config; kwargs...)
    flux = config.corrections isa SpinUp ? config.corrections.years : 0
    return greb_model!(RunSpec(flux=flux, ctrl=run.ctrl, scnr=run.scnr), _lower(config); kwargs...)
end
