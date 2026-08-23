"""
    PhysicsConfig

All model switches and parameters: mean-climate/CO₂-response toggles,
circulation components, hydrology parameters, external-forcing flags, and
the current experiment type. Passed explicitly to every physics function;
build one with [`create_experiment_config`](@ref) rather than the bare
keyword constructor for anything beyond `:full_model`.
"""
Base.@kwdef mutable struct PhysicsConfig
    # Mean Climate Switches
    log_clouds_dmc::Bool = true
    log_vapor_dmc::Bool = true
    log_crcl_dmc::Bool = true
    log_hydro_dmc::Bool = true
    log_atmos_dmc::Bool = true
    log_co2_dmc::Bool = true
    log_ocean_dmc::Bool = true
    log_qflux_dmc::Bool = true

    # CO₂ Response Switches
    log_clouds_drsp::Bool = true
    log_crcl_drsp::Bool = true
    log_hydro_drsp::Bool = true
    log_topo_drsp::Bool = true
    log_humid_drsp::Bool = true
    log_ocean_drsp::Bool = true

    # Circulation Components
    log_ice::Bool = true
    log_hdif::Bool = true
    log_hadv::Bool = true
    log_vdif::Bool = true
    log_vadv::Bool = true
    log_conv::Bool = true

    # Hydrology Parameters
    log_rain::Int = 0
    log_eva::Int = -1
    log_clim::Int = 0

    # External anomaly forcing gates
    log_tsurf_ext::Bool = false
    log_hwind_ext::Bool = false
    log_omega_ext::Bool = false

    # Experiment Type
    experiment::Symbol = :full_model  # :full_model, :constant_topo, :co2_double, etc.

    # CO₂ concentration for experiments (ppm)
    co2_concentration::Float32 = 340.0f0

    # Orbital-forcing experiments: which solar_scenarios table row to load
    # (:eccentricity / :obliquity, see load_solar_forcing_jld2)
    orbital_index::Int = 0

    # Earth-Sun distance experiment: percent change in orbital radius
    earth_sun_distance_pct::Float32 = 0.0f0

    # Hydrology parameters (calculated by set_hydrology_parameters!)
    c_q::Float32 = 1.0f0
    c_rq::Float32 = 0.0f0
    c_omega::Float32 = 0.0f0
    c_omegastd::Float32 = 0.0f0

    # IPCC scenario CO₂ lookup table (year => ppm)
    co2_scenario::Dict{Int,Float32} = Dict{Int,Float32}()

    # :custom_co2 experiment: path to a user-supplied "year CO2" text file
    # (same format as the IPCC scenario files), loaded into `co2_scenario`
    # at scenario start.
    custom_co2_path::String = ""
end

"""
    RunSpec

Run durations (in years) for [`greb_model!`](@ref): `flux` (flux-correction
spin-up), `ctrl` (control run), `scnr` (scenario run). A keyword struct
instead of three bare positional ints, whose order was easy to swap by
mistake.
"""
Base.@kwdef struct RunSpec
    flux::Int = 0
    ctrl::Int = 1
    scnr::Int = 1
end

# Static PhysicsConfig overrides per experiment. An empty NamedTuple means the
# experiment needs nothing beyond `experiment = sym`: `forcing` sets its CO₂
# per timestep, and `co2_concentration` is deliberately left at the default
# because it seeds the *control* run's CO₂ (see `init_model!`).
    :full_model             => (;),
    :constant_topo          => (log_topo_drsp = false,),
    :a1b_scenario           => (;),
    :co2_double             => (co2_concentration = 680.0f0,),
    :co2_quadruple          => (co2_concentration = 1360.0f0,),
    :co2_10x                => (;),
    :co2_half               => (;),
    :co2_zero               => (;),
    :co2_sine_wave          => (;),
    :co2_step               => (;),
    :solar_plus27           => (;),
    :solar_cycle_11yr       => (;),
    :paleo_231kyr           => (co2_concentration = 200.0f0,),
    :paleo_solar_modern_co2 => (;),
    :modern_solar_paleo_co2 => (;),
    :obliquity              => (;),
    :eccentricity           => (;),
    :earth_sun_distance     => (;),
    :elnino                 => (log_tsurf_ext = true, log_hwind_ext = true, log_omega_ext = true),
    :lanina                 => (log_tsurf_ext = true, log_hwind_ext = true, log_omega_ext = true),
    :rcp85                  => (log_tsurf_ext = true, log_hwind_ext = true, log_omega_ext = true),
    :rcp26                  => (;),
    :rcp45                  => (;),
    :rcp60                  => (;),
    :ssp119                 => (;),
    :ssp126                 => (;),
    :ssp245                 => (;),
    :ssp460                 => (;),
    :ssp585                 => (;),
    :historical_co2         => (;),
    :custom_co2             => (;),   # co2_path applied after construction
    :sst_plus1              => (;),
    :regional_co2_nh           => (;),
    :regional_co2_sh           => (;),
    :regional_co2_tropics      => (;),
    :regional_co2_extratropics => (;),
    :regional_co2_ocean        => (;),
    :regional_co2_land_ice     => (;),
    :regional_co2_winter       => (;),
    :regional_co2_summer       => (;),
)

# The 16 log_* keywords are consumed only by the two decon experiments.
const _DECON_ONLY_KEYWORDS = (:log_clouds_dmc, :log_ocean_dmc, :log_atmos_dmc,
    :log_co2_dmc, :log_hydro_dmc, :log_qflux_dmc, :log_topo_drsp,
    :log_clouds_drsp, :log_humid_drsp, :log_ocean_drsp, :log_hydro_drsp,
    :log_ice, :log_hdif, :log_hadv, :log_vdif, :log_vadv)

begin
    """
        create_experiment_config(experiment::Symbol; co2_path="", orbital_index=0,
            earth_sun_distance_pct=0.0, log_clouds_dmc=nothing, log_ocean_dmc=nothing,
            log_atmos_dmc=nothing, log_co2_dmc=nothing, log_hydro_dmc=nothing,
            log_qflux_dmc=nothing, log_topo_drsp=nothing, log_clouds_drsp=nothing,
            log_humid_drsp=nothing, log_ocean_drsp=nothing, log_hydro_drsp=nothing,
            log_ice=nothing, log_hdif=nothing, log_hadv=nothing, log_vdif=nothing,
            log_vadv=nothing) -> PhysicsConfig

    Create a pre-configured `PhysicsConfig` for any experiment the model
    dispatches on. Errors on an unknown symbol, listing the valid ones.

    The 16 `log_*` keywords are used **only** by
    `:decon_mean_climate`/`:decon_2xco2`; passing one for any other experiment
    warns and is ignored, rather than being silently dropped. `co2_path` applies
    only to `:custom_co2`, `orbital_index` to the paleo/orbital experiments, and
    `earth_sun_distance_pct` to `:earth_sun_distance`.

    Only `:constant_topo`, `:co2_double`, `:co2_quadruple`, `:paleo_231kyr` and
    the three forced-boundary experiments carry a static override; every other
    entry is `experiment = sym` alone, because [`forcing`](@ref) sets the
    scenario CO₂ per timestep and `co2_concentration` seeds the *control* run.

    # Experiments
    - `:full_model` - All processes active (default)
    - `:constant_topo` - Constant topography (log_topo_drsp = false), 2×CO₂ scenario
    - `:co2_double`/`:co2_quadruple`/`:co2_10x`/`:co2_half`/`:co2_zero` - CO₂ scaling
    - `:co2_sine_wave`/`:co2_step` - time-varying CO₂
    - `:a1b_scenario` - A1B CO₂ ramp (control baseline 280 ppm)
    - `:solar_plus27` - +27 W/m² solar constant
    - `:solar_cycle_11yr` - 11-year solar cycle
    - `:paleo_231kyr` - Paleoclimate (200 ppm CO₂)
    - `:paleo_solar_modern_co2`/`:modern_solar_paleo_co2` - crossed paleo/modern forcing
    - `:obliquity`/`:eccentricity` - orbital forcing, `orbital_index` selects the table row
    - `:earth_sun_distance` - solar constant scaled by `earth_sun_distance_pct`
    - `:elnino`/`:lanina` - ERA-Interim ENSO conditions
    - `:rcp26`/`:rcp45`/`:rcp60`/`:rcp85` - IPCC RCP climate change scenarios
    - `:ssp119`/`:ssp126`/`:ssp245`/`:ssp460`/`:ssp585` - IPCC SSP scenarios
    - `:historical_co2` - Observed CO₂ 1850-2017 (year starts at 1850, not 1950)
    - `:custom_co2` - user-supplied CO₂ trajectory, `co2_path` keyword gives
      the "year CO2" text file path (see [`load_custom_co2_scenario`](@ref))
    - `:sst_plus1` - ocean surface warmed 1 K, CO₂ held at control
    - `:regional_co2_nh`/`:regional_co2_sh`/`:regional_co2_tropics`/`:regional_co2_extratropics` -
      2×CO₂ over a latitude band (static mask, set by [`init_model!`](@ref))
    - `:regional_co2_ocean`/`:regional_co2_land_ice` - 2×CO₂ over ocean or
      land/ice (mask derived from the control run's ice cover, see
      [`apply_dynamic_co2_mask!`](@ref))
    - `:regional_co2_winter`/`:regional_co2_summer` - 2×CO₂ in one boreal season
    - `:decon_mean_climate` - deconstruct-mean-state experiment
    - `:decon_2xco2` - deconstruct-2×CO₂-response experiment
    """
    function create_experiment_config(experiment::Symbol; co2_path::AbstractString="",
        orbital_index::Int=0, earth_sun_distance_pct::Real=0.0,
        log_clouds_dmc::Union{Bool,Nothing}=nothing, log_ocean_dmc::Union{Bool,Nothing}=nothing,
        log_atmos_dmc::Union{Bool,Nothing}=nothing, log_co2_dmc::Union{Bool,Nothing}=nothing,
        log_hydro_dmc::Union{Bool,Nothing}=nothing, log_qflux_dmc::Union{Bool,Nothing}=nothing,
        log_topo_drsp::Union{Bool,Nothing}=nothing, log_clouds_drsp::Union{Bool,Nothing}=nothing,
        log_humid_drsp::Union{Bool,Nothing}=nothing, log_ocean_drsp::Union{Bool,Nothing}=nothing,
        log_hydro_drsp::Union{Bool,Nothing}=nothing, log_ice::Union{Bool,Nothing}=nothing,
        log_hdif::Union{Bool,Nothing}=nothing, log_hadv::Union{Bool,Nothing}=nothing,
        log_vdif::Union{Bool,Nothing}=nothing, log_vadv::Union{Bool,Nothing}=nothing)::PhysicsConfig

        on(x) = something(x, true)

        if experiment === :decon_mean_climate
            return PhysicsConfig(experiment=:decon_mean_climate,
                log_clouds_dmc=on(log_clouds_dmc), log_ocean_dmc=on(log_ocean_dmc),
                log_atmos_dmc=on(log_atmos_dmc), log_co2_dmc=on(log_co2_dmc),
                log_hydro_dmc=on(log_hydro_dmc), log_qflux_dmc=on(log_qflux_dmc),
                log_ice=on(log_ice), log_hdif=on(log_hdif), log_hadv=on(log_hadv),
                log_vdif=on(log_vdif), log_vadv=on(log_vadv))

        elseif experiment === :decon_2xco2
            return PhysicsConfig(experiment=:decon_2xco2, co2_concentration=680.0f0,
                log_topo_drsp=on(log_topo_drsp), log_clouds_drsp=on(log_clouds_drsp),
                log_humid_drsp=on(log_humid_drsp), log_ocean_drsp=on(log_ocean_drsp),
                log_hydro_drsp=on(log_hydro_drsp),
                log_ice=on(log_ice), log_hdif=on(log_hdif), log_hadv=on(log_hadv),
                log_vdif=on(log_vdif), log_vadv=on(log_vadv))
        end

        haskey(_EXPERIMENT_OVERRIDES, experiment) || error(
            "Unknown experiment: $experiment. Known: " *
            "$(sort!(vcat(collect(keys(_EXPERIMENT_OVERRIDES)), [:decon_mean_climate, :decon_2xco2])))")

        given = (log_clouds_dmc, log_ocean_dmc, log_atmos_dmc, log_co2_dmc,
            log_hydro_dmc, log_qflux_dmc, log_topo_drsp, log_clouds_drsp,
            log_humid_drsp, log_ocean_drsp, log_hydro_drsp, log_ice,
            log_hdif, log_hadv, log_vdif, log_vadv)
        passed = [kw for (kw, v) in zip(_DECON_ONLY_KEYWORDS, given) if v !== nothing]
        isempty(passed) || @warn "create_experiment_config(:$experiment) ignores " *
            "these keywords; only :decon_mean_climate/:decon_2xco2 use them. Set " *
            "the fields on the returned PhysicsConfig instead." keywords = passed

        cfg = PhysicsConfig(; experiment=experiment, _EXPERIMENT_OVERRIDES[experiment]...)
        experiment === :custom_co2 && (cfg.custom_co2_path = co2_path)
        experiment in (:obliquity, :eccentricity, :paleo_231kyr,
                       :paleo_solar_modern_co2) && (cfg.orbital_index = orbital_index)
        experiment === :earth_sun_distance &&
            (cfg.earth_sun_distance_pct = Float32(earth_sun_distance_pct))
        return cfg
    end
end;

# - Optimized Parameter Initialization ──────────────────────────────────
"""
    set_hydrology_parameters!(cfg::PhysicsConfig)

Initialize precipitation parameters `c_q, c_rq, c_omega, c_omegastd` based on
`cfg.log_rain` and `cfg.log_clim` settings.
"""
function set_hydrology_parameters!(cfg::PhysicsConfig)
    # Fast lookup instead of if-else chain
    haskey(HYDRO_PARAMS, cfg.log_rain) ||
        error("Unknown log_rain value: $(cfg.log_rain). Valid values: $(sort(collect(keys(HYDRO_PARAMS))))")
    params = HYDRO_PARAMS[cfg.log_rain]
    cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd = params

    # NCEP parameter adjustment
    if cfg.log_rain == 0 && cfg.log_clim == 1
        cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd = -1.27f0, 1.99f0, -16.54f0, 21.15f0
    end

    @info "⚙️ MSCM hydrology: log_rain=$(cfg.log_rain), log_clim=$(cfg.log_clim) → (c_q=$(cfg.c_q), c_rq=$(cfg.c_rq), c_omega=$(cfg.c_omega), c_omegastd=$(cfg.c_omegastd))"
end
