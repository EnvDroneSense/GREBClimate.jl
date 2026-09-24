"""
    apply_shortwave_addons!(state::ModelState, cfg, year, it) -> state.sw_solar_forcing

Multiplies the per-latitude shortwave multiplier `state.sw_solar_forcing` in
place by the add-ons set on `cfg`: the solar series `cfg.solar_scenario`, then
the aerosol `cfg.aerosol`. Unset add-ons are skipped, so default runs are
bit-identical. `year` and the scenario timestep `it` give the aerosol its
decimal year. Allocates nothing.
"""
function apply_shortwave_addons!(state::ModelState, cfg, year, it)
    mult = state.sw_solar_forcing
    if !isempty(cfg.solar_scenario)
        factor = get(cfg.solar_scenario, year, nothing)
        factor === nothing &&
            error("No solar data for year $year in cfg.solar_scenario " *
                  "(loaded $(length(cfg.solar_scenario)) years)")
        mult .*= factor
    end
    if cfg.aerosol !== nothing
        t = year + mod(it - 1, nstep_yr) / nstep_yr
        aerosol_transmission!(mult, state.aod, cfg.aerosol, t)
    end
    return mult
end
