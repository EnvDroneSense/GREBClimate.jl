# The CO2 and sunlight a scenario imposes at each timestep.

"""
    forcing(it, year, cfg::PhysicsConfig, fields::ClimateFields, icmn_ctrl; nstep_yr=nstep_yr)

Returns `(CO2, sw_solar_forcing)` for the current timestep, computed
according to `cfg.experiment`. Pure - the `regional_co2_*` masks are built
once per run by [`apply_dynamic_co2_mask!`](@ref), not here. `:full_model`
short-circuits before the experiment dispatch chain. The `:rcp26`/`:rcp45`/`:rcp60`/
`:custom_co2`/`:ssp*`/`:historical_co2` experiments look `year` up in
`cfg.co2_scenario`. With `cfg.log_co2_dmc` off, CO₂ is 0.
"""
function forcing(it, year, cfg::PhysicsConfig, fields::ClimateFields, icmn_ctrl; nstep_yr=nstep_yr)
    f = _experiment_forcing(it, year, cfg, nstep_yr)
    return cfg.log_co2_dmc ? f : (CO2=0.0f0, sw_solar_forcing=f.sw_solar_forcing)
end

function _experiment_forcing(it, year, cfg::PhysicsConfig, nstep_yr)
    # Default CO₂ concentration
    CO2 = cfg.co2_concentration
    sw_solar_forcing = 1.0f0

    # Fast path for the main experiment
    if cfg.experiment == :full_model
        return (CO2=CO2, sw_solar_forcing=sw_solar_forcing)
    end

    # - Legacy experiments ───────────
    if cfg.experiment == :constant_topo
        CO2 = 680.0f0  # 2x340 ppm

    elseif cfg.experiment == :a1b_scenario
        CO2_1950 = 310.0f0;
        CO2_2000 = 370.0f0;
        CO2_2050 = 520.0f0
        if year <= 2000
            CO2 = CO2_1950 + 60.0f0 / 50.0f0 * (year - 1950)
        elseif year <= 2050
            CO2 = CO2_2000 + 150.0f0 / 50.0f0 * (year - 2000)
        elseif year <= 2100
            CO2 = CO2_2050 + 180.0f0 / 50.0f0 * (year - 2050)
        end

    # - CO₂ scaling experiments ──────────────────────────────────────────────
    elseif cfg.experiment == :co2_double || cfg.experiment == :decon_2xco2
        CO2 = 680.0f0  # 2×CO₂

    elseif cfg.experiment == :co2_quadruple
        CO2 = 1360.0f0  # 4×CO₂

    elseif cfg.experiment == :co2_10x
        CO2 = 3400.0f0  # 10×CO₂

    elseif cfg.experiment == :co2_half
        CO2 = 170.0f0  # 0.5×CO₂

    elseif cfg.experiment == :co2_zero
        CO2 = 0.0f0  # 0×CO₂ (no greenhouse effect)

    # - Solar forcing experiments ───────────────────────────────────────────
    elseif cfg.experiment == :solar_plus27
        CO2 = 340.0f0
        sw_solar_forcing = (1365.0f0 + 27.0f0) / 1365.0f0

    elseif cfg.experiment == :solar_cycle_11yr
        CO2 = 340.0f0
        sw_solar_forcing = (1365.0f0 + 1.0f0 * sin(2f0*Float32(π) * year / 11.0f0)) / 1365.0f0

    # ── Time-varying CO₂ experiments ────────────
    elseif cfg.experiment == :co2_sine_wave
        CO2 = 340.0f0 + 170.0f0 + 170.0f0 * cos(2f0*Float32(π) * (year - 13.0f0) / 30.0f0)

    elseif cfg.experiment == :co2_step
        CO2 = year >= 1980 ? 340.0f0 : 680.0f0

    # ── Paleoclimate experiments ────────────────────
    elseif cfg.experiment == :paleo_231kyr
        CO2 = 200.0f0

    elseif cfg.experiment == :paleo_solar_modern_co2
        CO2 = 340.0f0

    elseif cfg.experiment == :modern_solar_paleo_co2
        CO2 = 200.0f0

    # ── Orbital forcing experiments ─────────────────
    elseif cfg.experiment == :obliquity
        CO2 = 340.0f0     # Solar forcing loaded externally

    elseif cfg.experiment == :eccentricity
        CO2 = 340.0f0     # Solar forcing loaded externally

    elseif cfg.experiment == :earth_sun_distance
        CO2 = 340.0f0     # Solar constant varies with Earth-Sun distance
        sw_solar_forcing = (1.0f0 / (1.0f0 + 0.01f0 * cfg.earth_sun_distance_pct))^2

    elseif cfg.experiment == :rcp85
        CO2 = 340.0f0  # Handled by boundary conditions

    # - IPCC RCP/SSP/historical/custom scenarios - CO₂ read from a per-year
    #   lookup table (`cfg.co2_scenario`, populated at scenario start) ───────
    elseif cfg.experiment in (:rcp26, :rcp45, :rcp60, :custom_co2,
                               :ssp119, :ssp126, :ssp245, :ssp460, :ssp585, :historical_co2)
        yr = round(Int, year)
        haskey(cfg.co2_scenario, yr) ||
            error("No CO2 data for year $yr in $(cfg.experiment) scenario table " *
                  "(loaded $(length(cfg.co2_scenario)) years)")
        CO2 = cfg.co2_scenario[yr]

    # - Regional/partial CO₂ experiments - static masks ─────────────────────
    elseif cfg.experiment in (:regional_co2_nh, :regional_co2_sh, :regional_co2_tropics, :regional_co2_extratropics)
        CO2 = 680.0f0

    # - Regional/partial CO₂ experiments - dynamic masks ────────────────────
    # `:regional_co2_ocean`/`:regional_co2_land_ice` need a mask derived from
    # the control run's ice cover; `apply_dynamic_co2_mask!` builds it once per
    # run in `greb_model!`, so nothing is computed per timestep here.
    elseif startswith(string(cfg.experiment), "regional_co2_")
        if cfg.experiment == :regional_co2_ocean
            # 2×CO₂ Ocean only - mask from apply_dynamic_co2_mask!
            CO2 = 680.0f0

        elseif cfg.experiment == :regional_co2_land_ice
            # 2×CO₂ Land/Ice only - mask from apply_dynamic_co2_mask!
            CO2 = 680.0f0

        elseif cfg.experiment == :regional_co2_winter
            # 2×CO₂ Boreal Winter only
            ityr_step = mod(it - 1, nstep_yr) + 1
            CO2 = (ityr_step <= 181 || ityr_step >= 547) ? 680.0f0 : 340.0f0

        elseif cfg.experiment == :regional_co2_summer
            # 2×CO₂ Boreal Summer only
            ityr_step = mod(it - 1, nstep_yr) + 1
            CO2 = (ityr_step <= 181 || ityr_step >= 547) ? 340.0f0 : 680.0f0
        end

        # - Forced boundary condition experiments (handled in scenario loop) ─────
    elseif cfg.experiment == :elnino || cfg.experiment == :lanina
        CO2 = 340.0f0
    end

    return (CO2=CO2, sw_solar_forcing=sw_solar_forcing)
end
