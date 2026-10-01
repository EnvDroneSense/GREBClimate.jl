# The net heat flux into the surface and into the atmosphere (W/m2): the sums
# both the flux-correction spin-up and the run integrate, defined once so a new
# term enters both.
#
# They are macros, not functions: a macro pastes the sum into each `@turbo`
# loop, where a function call changes how the loop contracts `em * LW_surf`
# into a fused multiply-add and so changes `Ta` in the last digit. The macros
# read the flux arrays by name from the calling scope: `SW`, `LW_surf`,
# `LW_down`, `Q_lat`, `Q_sens` for the surface; `LW_up`, `LW_down`, `em`,
# `LW_surf`, `Q_lat_air`, `Q_sens` for the atmosphere.

macro surface_flux(i, j)
    return esc(:(SW[$i, $j] + LW_surf[$i, $j] - LW_down[$i, $j] + Q_lat[$i, $j] + Q_sens[$i, $j]))
end

macro atmosphere_flux(i, j)
    return esc(:(LW_up[$i, $j] + LW_down[$i, $j] - em[$i, $j] * LW_surf[$i, $j] + Q_lat_air[$i, $j] - Q_sens[$i, $j]))
end

# The same sums at one cell of a `tendencies!` result, for tests and tools.
# Outside `@turbo` they can differ from the loops' values in the last digit.
function surface_flux(tend, i, j)
    SW, LW_surf, LW_down, Q_lat, Q_sens = tend.SW, tend.LW_surf, tend.LW_down, tend.Q_lat, tend.Q_sens
    return @surface_flux(i, j)
end

function atmosphere_flux(tend, i, j)
    LW_up, LW_down, em, LW_surf, Q_lat_air, Q_sens = tend.LW_up, tend.LW_down, tend.em, tend.LW_surf, tend.Q_lat_air, tend.Q_sens
    return @atmosphere_flux(i, j)
end
