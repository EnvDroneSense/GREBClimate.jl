# The CO2 and sunlight a scenario imposes at each timestep, and where its CO2
# applies.

"""
    forcing(t::ModelTime, r::ResolvedConfig) -> (CO2, sw_solar_forcing)

The scenario's CO2 (ppm) and solar multiplier at time `t`, from its
[`CO2Path`](@ref) and [`Solar`](@ref) parts. CO2 is
0 with `Processes(co2 = false)`. Pure: where the CO2 applies is set once, at
the start of the scenario, by `apply_co2_mask!` and `apply_surface_mask!`.
"""
function forcing(t::ModelTime, r::ResolvedConfig)
    s = r.config.scenario
    co2 = r.config.processes.co2 ? _co2_at(s.co2, t, r) : 0.0f0
    return (CO2=co2, sw_solar_forcing=_solar_factor(s.solar, t))
end

_co2_at(c::ConstantCO2, t, r) = c.ppm

function _co2_at(::Union{CO2Table,CO2File}, t, r)
    haskey(r.co2_table, year(t)) ||
        error("No CO2 data for year $(year(t)) in the scenario's CO2 table (loaded $(length(r.co2_table)) years)")
    return r.co2_table[year(t)]
end

# After 2100 the ramp holds its last value
function _co2_at(::A1BRamp, t, r)
    y = year(t)
    CO2_1950 = 310.0f0
    CO2_2000 = 370.0f0
    CO2_2050 = 520.0f0
    if y <= 2000
        return CO2_1950 + 60.0f0 / 50.0f0 * (y - 1950)
    elseif y <= 2050
        return CO2_2000 + 150.0f0 / 50.0f0 * (y - 2000)
    elseif y <= 2100
        return CO2_2050 + 180.0f0 / 50.0f0 * (y - 2050)
    end
    return 700.0f0
end

_co2_at(::CO2SineWave, t, r) = 510.0f0 + 170.0f0 * cos(2f0*Float32(π) * (year(t) - 13.0f0) / 30.0f0)

_co2_at(c::CO2Step, t, r) = year(t) >= c.year ? c.after : c.before

# Boreal winter is half the year: from the first step of 1 October through the
# first step of 1 April, as in the original code
const _winter_first_step = first_step_of(10, 1)
const _winter_last_step = first_step_of(4, 1)

function _co2_at(c::SeasonalCO2, t, r)
    step = step_of_year(t)
    winter = step <= _winter_last_step || step >= _winter_first_step
    return winter == (c.season === :boreal_winter) ? c.in_season : c.out_of_season
end

_solar_factor(::Union{ModernSolar,SolarTable}, t) = 1.0f0
_solar_factor(s::SolarConstant, t) = (1365.0f0 + s.offset) / 1365.0f0
_solar_factor(s::SolarCycle, t) = (1365.0f0 + s.amplitude * sin(2f0*Float32(π) * year(t) / s.period)) / 1365.0f0
_solar_factor(s::EarthSunDistance, t) = (1.0f0 / (1.0f0 + 0.01f0 * s.percent))^2

"""
    apply_co2_mask!(mask::CO2Mask, fields::ClimateFields)

Sets `fields.co2_part`, the fraction of the scenario CO2 each cell gets: 1
everywhere, then 0.5 outside a [`LatitudeMask`](@ref) band. `greb_model!`
calls it at the start of the scenario; the spin-up and the control run on the
full CO2 everywhere. A [`SurfaceMask`](@ref) needs the control run's ice cover
and is set by `apply_surface_mask!`.
"""
function apply_co2_mask!(mask::CO2Mask, fields::ClimateFields)
    fields.co2_part .= 1.0f0
    _latitude_mask!(fields.co2_part, mask)
    return nothing
end

_latitude_mask!(co2_part, ::CO2Mask) = co2_part

# The tropics band runs from 33.75 S to 30 N, as in the original code. In the
# two rows at the band edge that are halved, every fourth longitude keeps the
# full CO2.
const _tropics_south = -33.75f0
const _tropics_north = 30.0f0

function _latitude_mask!(co2_part, m::LatitudeMask)
    south = findall(<(0), lat_grid)
    tropics = findall(lat -> _tropics_south < lat < _tropics_north, lat_grid)
    if m.band === :nh
        co2_part[:, south] .= 0.5f0
    elseif m.band === :sh
        co2_part[:, setdiff(1:ydim, south)] .= 0.5f0
    elseif m.band === :tropics
        co2_part[:, setdiff(1:ydim, tropics)] .= 0.5f0
        co2_part[4:4:xdim, [first(tropics) - 1, last(tropics) + 1]] .= 1.0f0
    else
        co2_part[:, tropics] .= 0.5f0
        co2_part[4:4:xdim, [first(tropics), last(tropics)]] .= 1.0f0
    end
    return co2_part
end

"""
    apply_surface_mask!(mask::CO2Mask, fields::ClimateFields, icmn_ctrl)

Sets `fields.co2_part` for a [`SurfaceMask`](@ref) from the control run's
annual-mean ice cover `icmn_ctrl`: `:ocean` halves CO2 over land and over ice,
`:land_ice` halves it over ice-free ocean. A no-op for every other mask.
"""
apply_surface_mask!(::CO2Mask, fields::ClimateFields, icmn_ctrl) = nothing

function apply_surface_mask!(mask::SurfaceMask, fields::ClimateFields, icmn_ctrl)
    co2_part = fields.co2_part
    z_topo = fields.z_topo
    co2_part .= 1.0f0

    # Annual-mean ice cover, not month 1
    icmn_ctrl1 = dropdims(sum(icmn_ctrl, dims=3), dims=3) ./ size(icmn_ctrl, 3)

    if mask.surface === :ocean
        # 2×CO2 ocean only: halve CO2 over land, and over annual-mean ice.
        for j in 1:ydim, i in 1:xdim
            if is_land(z_topo[i, j])
                co2_part[i, j] = 0.5f0
            end
        end
        for j in 1:ydim, i in 1:xdim
            if icmn_ctrl1[i, j] >= 0.5f0
                co2_part[i, j] = 0.5f0
            end
        end
    else
        # 2×CO2 land/ice only: halve CO2 over ocean, then exempt annual-mean ice.
        for j in 1:ydim, i in 1:xdim
            if !is_land(z_topo[i, j])
                co2_part[i, j] = 0.5f0
            end
        end
        for j in 1:ydim, i in 1:xdim
            if icmn_ctrl1[i, j] >= 0.5f0
                co2_part[i, j] = 1.0f0
            end
        end
    end
    return nothing
end
