# How the annual summary line names the year: the spin-up counts its years from 1
_year_label(t::ModelTime) =
    t.phase === spinup ? "spin-up year $(year(t) - t.start_year + 1)" : string(year(t))

# The two cells of the annual summary line: label, longitude index, latitude index
const _SAMPLE_CELLS = (("178 E 9 N", 48, 27), ("58 E 51 N", 16, 38))

"""
    diagnostics!(t::ModelTime, surf::SurfaceState, state)

Accumulates the current timestep into `state`'s annual-mean buffers; at the
last timestep of the year, averages them, logs the annual summary line
(the area-weighted global mean and two sample cells, in °C) with `@info`, and
resets the accumulators for the next year.
"""
function diagnostics!(t::ModelTime, surf::SurfaceState, state::ModelState)
    # Accumulate
    Ts_annual_mean = state.Ts_annual_mean
    Ts = surf.Ts

    @turbo for j in 1:ydim
        for i in 1:xdim
            Ts_annual_mean[i, j] += Ts[i, j]
        end
    end

    if is_year_end(t)
        # Compute annual means
        n = nstep_yr
        state.Ts_annual_mean ./= n

        # Annual-mean surface temperature (°C): global mean and the two sample cells
        celsius(T) = round(T - 273.15; digits=2)
        cells = join(("$label $(celsius(state.Ts_annual_mean[i, j]))" for (label, i, j) in _SAMPLE_CELLS), "; ")
        @info "$(_year_label(t)): Ts global mean $(celsius(global_mean(state.Ts_annual_mean))) °C; $cells"

        # Reset accumulators
        fill!(state.Ts_annual_mean, 0.0f0)
    end
    return nothing
end

"""
    output!(t::ModelTime, surf::SurfaceState, tend, ws, output_buf, times, acc)

Accumulates the current timestep into `acc`; on the last timestep of a month,
pushes the monthly-mean [`MonthlyRecord`](@ref) onto `output_buf` and its
[`RecordTime`](@ref) onto `times`, and resets `acc`. `tend` is the
`NamedTuple` [`tendencies!`](@ref) returns; `ws.precip`/`evap`/
`qcrcl` hold this step's converted precipitation/evaporation/moisture-
circulation output.
"""
function output!(t::ModelTime, surf::SurfaceState, tend, ws::ModelWorkspace,
    output_buf::Vector{MonthlyRecord}, times::Vector{RecordTime}, acc::MonthlyAccumulator)
    accumulate!(acc, surf, tend, ws)

    if is_month_end(t)
        push!(output_buf, monthly_means(acc, steps_in_month(month(t))))
        push!(times, (year=year(t), month=month(t)))
        reset!(acc)
    end
    return nothing
end

"""
    time_loop!(t::ModelTime, CO2, Ts, Ta, q, To, output_buf, times, fields, state, ws, acc, timestate, r::ResolvedConfig;
               ws_a=ws, ws_q=ws, observer=nothing)

One full model timestep at time `t`: computes [`tendencies!`](@ref),
integrates `Ts`/`Ta`/`To`/`q` forward with flux corrections applied, runs
[`seaice!`](@ref), then dispatches to [`output!`](@ref) and
[`diagnostics!`](@ref). `ws_a`/`ws_q` are forwarded to
[`tendencies!`](@ref) - see its docstring for the opt-in threading they
enable. `observer` is called before and after the update (see
[`greb_model!`](@ref)).
"""
function time_loop!(t::ModelTime, CO2, Ts, Ta, q, To, output_buf, times,
    fields::ClimateFields, state::ModelState, ws::ModelWorkspace, acc::MonthlyAccumulator,
    timestate, r::ResolvedConfig; ws_a::ModelWorkspace=ws, ws_q::ModelWorkspace=ws,
    observer=nothing)
    timestate.jday = day_of_year(t)
    timestate.ityr = step_of_year(t)
    ityr = timestate.ityr

    # Compute tendencies
    tend = tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, timestate, r; ws_a=ws_a, ws_q=ws_q)

    observer === nothing ||
        observer(:after_tendencies, _step_view(t, CO2, Ts, Ta, To, q, tend, fields, r))

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
        observer(:after_step, _step_view(t, CO2, Ts, Ta, To, q, tend, fields, r))

    # Output and diagnostics
    surf = SurfaceState(Ts, Ta, To, q)
    output!(t, surf, tend, ws, output_buf, times, acc)
    diagnostics!(t, surf, state)

    return nothing
end
