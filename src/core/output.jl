# What the annual summary line of a spin-up year shows in place of a calendar year
struct _SpinUpYear
    n::Int
end
Base.show(io::IO, y::_SpinUpYear) = print(io, "spin-up year ", y.n)

# The two cells of the annual summary line: label, longitude index, latitude index
const _SAMPLE_CELLS = (("178 E 9 N", 48, 27), ("58 E 51 N", 16, 38))

"""
    diagnostics!(year, surf::SurfaceState, state, timestate)

Accumulates the current timestep into `state`'s annual-mean buffers; at the
last timestep of the year, averages them, logs the annual summary line
(the area-weighted global mean and two sample cells, in °C) with `@info`, and
resets the accumulators for the next year.
"""
function diagnostics!(year, surf::SurfaceState, state::ModelState, timestate)
    # Accumulate
    Ts_annual_mean = state.Ts_annual_mean
    Ts = surf.Ts

    @turbo for j in 1:ydim
        for i in 1:xdim
            Ts_annual_mean[i, j] += Ts[i, j]
        end
    end

    if timestate.ityr == nstep_yr
        # Compute annual means
        n = nstep_yr
        state.Ts_annual_mean ./= n

        # Annual-mean surface temperature (°C): global mean and the two sample cells
        celsius(T) = round(T - 273.15; digits=2)
        cells = join(("$label $(celsius(state.Ts_annual_mean[i, j]))" for (label, i, j) in _SAMPLE_CELLS), "; ")
        @info "$year: Ts global mean $(celsius(global_mean(state.Ts_annual_mean))) °C; $cells"

        # Reset accumulators
        fill!(state.Ts_annual_mean, 0.0f0)
    end
    return nothing
end

"""
    output!(it, irec, mon, surf::SurfaceState, tend, ws, output_buf, acc, timestate)

Accumulates the current timestep into `acc`; on the last timestep of `mon`,
pushes a monthly-mean [`MonthlyRecord`](@ref) onto `output_buf`, resets `acc`,
and advances to the next month. Returns `(mon, irec)`. `tend` is the
`NamedTuple` [`tendencies!`](@ref) returns; `ws.precip`/`evap`/
`qcrcl` hold this step's converted precipitation/evaporation/moisture-
circulation output.
"""
function output!(it, irec, mon, surf::SurfaceState, tend, ws::ModelWorkspace,
    output_buf::Vector{MonthlyRecord}, acc::MonthlyAccumulator, timestate)
    mon = clamp(mon, 1, months_per_year)

    accumulate!(acc, surf.Ts, surf.Ta, surf.To, surf.q, tend.albedo, tend.ice_cover,
        ws.precip, ws.evap, ws.qcrcl, tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens,
        tend.LW_up, tend.LW_down, tend.em)

    # ----- Check end of month -----
    if timestate.jday == jday_mon_cumsum[mon] && is_day_end(it)
        ndm = steps_in_month(mon)
        irec += 1
        push!(output_buf, (
            Ts=acc.Tmm ./ ndm,
            Ta=acc.Tamm ./ ndm,
            To=acc.Tomm ./ ndm,
            q=acc.qmm ./ ndm,
            albedo=acc.apmm ./ ndm,
            ice=acc.icemm ./ ndm,
            precip=acc.precipmm ./ ndm,
            evap=acc.evapmm ./ ndm,
            qcrcl=acc.qcrclmm ./ ndm,
            sw=acc.swmm ./ ndm,
            lw=acc.lwmm ./ ndm,
            qlat=acc.qlatmm ./ ndm,
            qsens=acc.qsensmm ./ ndm,
            olr=acc.olrmm ./ ndm,
            lwdown=acc.lwdownmm ./ ndm
        ))
        reset!(acc)
        mon = mod(mon, months_per_year) + 1
    end
    return (mon=mon, irec=irec)
end

"""
    time_loop!(it, year, CO2, mon, irec, Ts, Ta, q, To, output_buf, fields, state, ws, acc, timestate, r::ResolvedConfig;
               ws_a=ws, ws_q=ws, observer=nothing, phase=:ctrl)

One full model timestep: computes [`tendencies!`](@ref), integrates
`Ts`/`Ta`/`To`/`q` forward with flux corrections applied, runs
[`seaice!`](@ref), then dispatches to [`output!`](@ref) and
[`diagnostics!`](@ref). Returns `(mon, irec)`. `ws_a`/`ws_q` are forwarded to
[`tendencies!`](@ref) - see its docstring for the opt-in threading they
enable. `observer` is called before and after the update (see
[`greb_model!`](@ref)); `phase` is passed on to it.
"""
function time_loop!(it, year, CO2, mon, irec, Ts, Ta, q, To, output_buf,
    fields::ClimateFields, state::ModelState, ws::ModelWorkspace, acc::MonthlyAccumulator,
    timestate, r::ResolvedConfig; ws_a::ModelWorkspace=ws, ws_q::ModelWorkspace=ws,
    observer=nothing, phase::Symbol=:ctrl)
    timestate.jday = day_of_year(it)
    timestate.ityr = step_of_year(it)
    ityr = timestate.ityr

    # Compute tendencies
    tend = tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, timestate, r; ws_a=ws_a, ws_q=ws_q)

    observer === nothing ||
        observer(:after_tendencies, _step_view(phase, it, year, ityr, CO2, Ts, Ta, To, q, tend, fields, r))

    # Correction views
    TF_corr = @view fields.Ts_flux_correction[:, :, ityr]
    qF_corr = @view fields.q_flux_correction[:, :, ityr]
    ToF_corr = @view fields.To_flux_correction[:, :, ityr]
    cap_surf = fields.cap_surf
    wz_vapor = fields.wz_vapor

    # Humidity tendency buffer selection
    dq_eva_use = tend.dq_eva
    dq_rain_use = tend.dq_rain
    dq_crcl_use = tend.dq_crcl
    hydro_on = r.config.processes.hydrology !== :none ? 1.0f0 : 0.0f0

    SW = tend.SW; LW_surf = tend.LW_surf; LW_down = tend.LW_down
    Q_lat = tend.Q_lat; Q_sens = tend.Q_sens; dTa_crcl = tend.dTa_crcl
    LW_up = tend.LW_up; em = tend.em; Q_lat_air = tend.Q_lat_air
    dTo = tend.dTo; dT_ocean = tend.dT_ocean
    precip = ws.precip; evap = ws.evap; qcrcl = ws.qcrcl

    # Surface/air temperature, deep ocean, and humidity update
    @turbo for j in 1:ydim
        for i in 1:xdim
            Ts[i, j] = Ts[i, j] + dT_ocean[i, j] + Δt * (@surface_flux(i, j) + TF_corr[i, j]) / cap_surf[i, j]
            Ta[i, j] = Ta[i, j] + dTa_crcl[i, j] + Δt * @atmosphere_flux(i, j) / cap_air

            Ts[i, j] = ifelse(Ts[i, j] < min_T_K, min_T_K, Ts[i, j])
            Ta[i, j] = ifelse(Ta[i, j] < min_T_K, min_T_K, Ta[i, j])

            To[i, j] = To[i, j] + dTo[i, j] + ToF_corr[i, j]

            tb = Δt * (dq_eva_use[i, j] + dq_rain_use[i, j]) + dq_crcl_use[i, j] + qF_corr[i, j]
            tb = ifelse(tb <= -q[i, j], -min_humidity_change * q[i, j], tb)
            tb = ifelse(tb > max_humidity_change, max_humidity_change, tb)
            tb = hydro_on * tb
            q[i, j] = q[i, j] + tb

            precip[i, j] = (-dq_rain_use[i, j]) * wz_vapor[i, j] * q_to_mm_per_day
            evap[i, j] = dq_eva_use[i, j] * wz_vapor[i, j] * q_to_mm_per_day
            qcrcl[i, j] = dq_crcl_use[i, j]
        end
    end

    # Sea ice heat capacity
    seaice!(Ts, fields, timestate, r.config.processes)

    observer === nothing ||
        observer(:after_step, _step_view(phase, it, year, ityr, CO2, Ts, Ta, To, q, tend, fields, r))

    # Output and diagnostics
    surf = SurfaceState(Ts, Ta, To, q)
    (; mon, irec) = output!(it, irec, mon, surf, tend, ws, output_buf, acc, timestate)
    diagnostics!(year, surf, state, timestate)

    return (mon=mon, irec=irec)
end
