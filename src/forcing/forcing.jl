# The CO2 and sunlight a scenario imposes at each timestep, and where its CO2
# applies.

"""
    forcing(it, year, r::ResolvedConfig) -> (CO2, sw_solar_forcing)

The scenario's CO₂ (ppm) and solar multiplier at scenario step `it` in
calendar `year`, from its [`CO2Path`](@ref) and [`Solar`](@ref) parts. CO₂ is
0 with `Processes(co2 = false)`. Pure: where the CO₂ applies is set once per
run by `apply_co2_mask!` and `apply_dynamic_co2_mask!`.
"""
function forcing(it, year, r::ResolvedConfig)
    s = r.config.scenario
    co2 = r.config.processes.co2 ? _co2_at(s.co2, it, year, r) : 0.0f0
    return (CO2=co2, sw_solar_forcing=_solar_factor(s.solar, year))
end

_co2_at(c::ConstantCO2, it, year, r) = c.ppm

function _co2_at(::Union{CO2Table,CO2File}, it, year, r)
    yr = round(Int, year)
    haskey(r.co2_table, yr) ||
        error("No CO2 data for year $yr in the scenario's CO2 table (loaded $(length(r.co2_table)) years)")
    return r.co2_table[yr]
end

# After 2100 the ramp falls back to 340 ppm, as the original code does
function _co2_at(::A1BRamp, it, year, r)
    CO2_1950 = 310.0f0
    CO2_2000 = 370.0f0
    CO2_2050 = 520.0f0
    if year <= 2000
        return CO2_1950 + 60.0f0 / 50.0f0 * (year - 1950)
    elseif year <= 2050
        return CO2_2000 + 150.0f0 / 50.0f0 * (year - 2000)
    elseif year <= 2100
        return CO2_2050 + 180.0f0 / 50.0f0 * (year - 2050)
    end
    return 340.0f0
end

_co2_at(::CO2SineWave, it, year, r) = 340.0f0 + 170.0f0 + 170.0f0 * cos(2f0*Float32(π) * (year - 13.0f0) / 30.0f0)

_co2_at(c::CO2Step, it, year, r) = year >= c.year ? c.after : c.before

function _co2_at(c::SeasonalCO2, it, year, r)
    step = mod(it - 1, nstep_yr) + 1
    winter = step <= 181 || step >= 547
    return winter == (c.season === :boreal_winter) ? c.inside : c.outside
end

_solar_factor(::Union{ModernSolar,SolarTable}, year) = 1.0f0
_solar_factor(s::SolarConstant, year) = (1365.0f0 + s.dW) / 1365.0f0
_solar_factor(s::SolarCycle, year) = (1365.0f0 + s.amplitude * sin(2f0*Float32(π) * year / s.period)) / 1365.0f0
_solar_factor(s::EarthSunDistance, year) = (1.0f0 / (1.0f0 + 0.01f0 * s.pct))^2

"""
    apply_co2_mask!(mask::CO2Mask, fields::ClimateFields)

Sets `fields.co2_part`, the fraction of the scenario CO₂ each cell gets: 1
everywhere, then 0.5 outside a [`LatitudeMask`](@ref) band. A
[`SurfaceMask`](@ref) needs the control run's ice cover and is set later by
`apply_dynamic_co2_mask!`.
"""
function apply_co2_mask!(mask::CO2Mask, fields::ClimateFields)
    fields.co2_part .= 1.0f0
    _latitude_mask!(fields.co2_part, mask)
    return nothing
end

_latitude_mask!(co2_part, ::CO2Mask) = co2_part

# The band edges at every fourth longitude follow the original code
function _latitude_mask!(co2_part, m::LatitudeMask)
    if m.band === :nh
        co2_part[:, 1:24] .= 0.5f0
    elseif m.band === :sh
        co2_part[:, 25:48] .= 0.5f0
    elseif m.band === :tropics
        co2_part[:, 1:15] .= 0.5f0
        co2_part[:, 33:48] .= 0.5f0
        for i in 4:4:96
            co2_part[i, 33] = 1.0f0
            co2_part[i, 15] = 1.0f0
        end
    else
        co2_part[:, 16:32] .= 0.5f0
        for i in 4:4:96
            co2_part[i, 32] = 1.0f0
            co2_part[i, 16] = 1.0f0
        end
    end
    return co2_part
end

"""
    apply_dynamic_co2_mask!(mask::CO2Mask, fields::ClimateFields, icmn_ctrl)

Sets `fields.co2_part` for a [`SurfaceMask`](@ref) from the control run's
annual-mean ice cover `icmn_ctrl`: `:ocean` halves CO₂ over land and over ice,
`:land_ice` halves it over ice-free ocean. A no-op for every other mask.
"""
apply_dynamic_co2_mask!(::CO2Mask, fields::ClimateFields, icmn_ctrl) = nothing

function apply_dynamic_co2_mask!(mask::SurfaceMask, fields::ClimateFields, icmn_ctrl)
    co2_part = fields.co2_part
    z_topo = fields.z_topo
    co2_part .= 1.0f0

    # Annual-mean ice cover, not month 1
    icmn_ctrl1 = dropdims(sum(icmn_ctrl, dims=3), dims=3) ./ size(icmn_ctrl, 3)

    if mask.surface === :ocean
        # 2×CO₂ ocean only: halve CO₂ over land, and over annual-mean ice.
        for j in 1:ydim, i in 1:xdim
            if z_topo[i, j] > 0.0f0
                co2_part[i, j] = 0.5f0
            end
        end
        for j in 1:ydim, i in 1:xdim
            if icmn_ctrl1[i, j] >= 0.5f0
                co2_part[i, j] = 0.5f0
            end
        end
    else
        # 2×CO₂ land/ice only: halve CO₂ over ocean, then exempt annual-mean ice.
        for j in 1:ydim, i in 1:xdim
            if z_topo[i, j] <= 0.0f0
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
