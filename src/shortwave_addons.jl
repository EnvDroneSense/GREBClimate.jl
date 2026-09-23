"""
    apply_shortwave_addons!(mult, cfg, year, it) -> mult

Multiply the per-latitude shortwave multiplier `mult` in place by the optional
add-ons configured on `cfg`. Currently applies `cfg.solar_scenario` (a year to
multiplier table); a no-op when it is empty, so default runs are bit-identical.
`it` is the timestep index, unused here but reserved for add-ons that need
sub-year time.
"""
function apply_shortwave_addons!(mult::Vector{Float32}, cfg, year, it)
    if !isempty(cfg.solar_scenario)
        factor = get(cfg.solar_scenario, year, nothing)
        factor === nothing &&
            error("No solar data for year $year in cfg.solar_scenario " *
                  "(loaded $(length(cfg.solar_scenario)) years)")
        mult .*= factor
    end
    return mult
end
