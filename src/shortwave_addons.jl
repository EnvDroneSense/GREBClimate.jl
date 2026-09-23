"""
    apply_shortwave_addons!(state::ModelState, cfg, year, it) -> state.sw_solar_forcing

Multiply the per-latitude shortwave multiplier `state.sw_solar_forcing` in place
by the optional add-ons configured on `cfg`, in this order: the solar series
`cfg.solar_scenario` (a year to multiplier table), then the stratospheric aerosol
`cfg.aerosol`. Each is skipped when unset, so default runs are bit-identical.
`year` is the integer scenario year and `it` the scenario timestep, which
together give the aerosol its sub-year time. Allocates nothing.
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
