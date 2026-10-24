# What the `Processes` options and the correction mode replace in `fields`: one
# row per override, `(applies, field, set!)`, grouped by field and applied in
# order. `init_model!` applies the rows that apply to a run and `_mutated_fields`
# names their fields, so the two cannot disagree
const _OVERRIDES = (
    (c -> c.processes.topography === :flat, :z_topo, a -> (@. a = min(a, 1.0f0))),   # constant topography
    (c -> c.processes.clouds === :none, :cloud_clim, a -> fill!(a, 0.0f0)),
    (c -> c.processes.clouds === :uniform, :cloud_clim, a -> fill!(a, 0.7f0)),       # 2xCO2 deconstruction
    (c -> c.corrections isa NoCorrections, :Ts_flux_correction, a -> fill!(a, 0.0f0)),
    (c -> c.corrections isa NoCorrections, :q_flux_correction, a -> fill!(a, 0.0f0)),
    (c -> c.corrections isa NoCorrections, :To_flux_correction, a -> fill!(a, 0.0f0)),
    (c -> c.processes.hydrology === :none, :q_clim, a -> fill!(a, 0.0f0)),
    (c -> c.processes.humidity === :uniform, :q_clim, a -> fill!(a, 0.0052f0)),      # constant water vapor
    (c -> c.processes.ocean === :mixed_layer, :mld_clim, a -> fill!(a, d_ocean)),    # no deep ocean
)

"""
    init_model!(r::ResolvedConfig, fields::ClimateFields)

One-time per-run setup: resets the CO2 mask to full CO2, applies the climatologies the
[`Processes`](@ref) options replace (flat topography, cloud cover, humidity,
mixed layer) and zeroes the flux corrections for `NoCorrections()`, derives
the fields that follow from them ([`derive_fields!`](@ref)), then computes the
control-run initial state. Returns `(surf, CO2_ctrl)`: `surf` is the
[`SurfaceState`](@ref) at the last step of the climatology (the air starts at the
surface temperature), `CO2_ctrl` the control CO2.
"""
function init_model!(r::ResolvedConfig, fields::ClimateFields)
    c = r.config
    p = c.processes

    # The spin-up and the control run on the full CO2 everywhere; a regional
    # mask is set at the start of the scenario. The reset also clears what an
    # earlier regional run left in a reused `fields`
    fields.co2_part .= 1.0f0

    for (applies, field, set!) in _OVERRIDES
        applies(c) && set!(getfield(fields, field))
    end

    derive_fields!(fields, p)

    # Initial conditions from the last time step of the climatology
    Ts = fields.Ts_clim[:, :, nstep_yr] |> copy
    surf = SurfaceState(Ts, copy(Ts), fields.To_clim[:, :, nstep_yr] |> copy,
        fields.q_clim[:, :, nstep_yr] |> copy)

    CO2_ctrl = p.co2 ? c.scenario.control_co2 : 0.0f0

    return (surf=surf, CO2_ctrl=CO2_ctrl)
end

_add_boundary_anomaly!(::SurfaceForcing, fields::ClimateFields) = fields

function _add_boundary_anomaly!(b::BoundaryAnomaly, fields::ClimateFields)
    anomaly = fields.boundary_anomaly
    (anomaly === nothing || fields.boundary_anomaly_source[2] !== b.source) &&
        error("the boundary anomaly :$(b.source) is not loaded: call load_boundary_anomaly! first")
    if b.source === :cmip5_rcp85
        @info "Applying CMIP5 RCP8.5 climate change forcing"
    else
        @info "Applying ERA-Interim $(b.source === :elnino ? "El Niño" : "La Niña") forcing"
    end
    for name in _BOUNDARY_FIELDS
        getfield(fields, name) .+= anomaly[name]
    end
    return fields
end

"""
    _mutated_fields(c::Config) -> Vector{Symbol}

The arrays of `fields` a run of `c` overwrites and [`greb_model!`](@ref)
restores: the solar table, the flux corrections, the boundary climatologies
(for a [`BoundaryAnomaly`](@ref)) and the climatologies the
[`Processes`](@ref) options replace. Everything else `init_model!` writes is
derived from these.
"""
function _mutated_fields(c::Config)
    names = [:sw_solar, :Ts_flux_correction, :q_flux_correction, :To_flux_correction]
    c.scenario.surface isa BoundaryAnomaly && append!(names, _BOUNDARY_FIELDS)
    for (applies, field, _) in _OVERRIDES
        applies(c) && !(field in names) && push!(names, field)
    end
    return names
end

# The control run and the spin-up carry this year as their label
const _control_start_year = 1970

"""
    qflux_correction!(CO2_ctrl, Ts, Ta, To, q, fields, state, r::ResolvedConfig, ws, years; ws_a=ws, ws_q=ws)

Runs `years` years of `tendencies!` to derive the ocean/atmosphere flux
corrections (`fields.Ts_flux_correction`/`q_flux_correction`/`To_flux_correction`) that make the
control climate match observed climatology. Mutates `Ts`/`Ta`/`To`/`q` in
place as it integrates. `ws_a`/`ws_q` are forwarded to [`tendencies!`](@ref).
"""
function qflux_correction!(CO2_ctrl, Ts, Ta, To, q, fields::ClimateFields, state::ModelState, r::ResolvedConfig, ws::ModelWorkspace, years;
    ws_a::ModelWorkspace=ws, ws_q::ModelWorkspace=ws)
    cap_surf = fields.cap_surf
    for t in steps(spinup, _control_start_year, years)
        slice = data_slice(t)

        tend = tendencies!(CO2_ctrl, Ts, Ta, To, q, fields, state, ws, slice, r; ws_a=ws_a, ws_q=ws_q)

        Tc = clim_slice(fields.Ts_clim, slice)
        Toc = clim_slice(fields.To_clim, slice)
        qc = clim_slice(fields.q_clim, slice)
        TFc = clim_slice(fields.Ts_flux_correction, slice)
        ToFc = clim_slice(fields.To_flux_correction, slice)
        qFc = clim_slice(fields.q_flux_correction, slice)

        # Surface/air temperature, deep ocean, and humidity update.
        Ts0_buf = ws.Ts0; Ta0_buf = ws.Ta0; To0_buf = ws.To0; q0_buf = ws.q0
        dT_ocean = tend.dT_ocean; SW = tend.SW; LW_surf = tend.LW_surf; LW_down = tend.LW_down
        Q_lat = tend.Q_lat; Q_sens = tend.Q_sens; dTa_crcl = tend.dTa_crcl
        LW_up = tend.LW_up; em = tend.em; Q_lat_air = tend.Q_lat_air
        dTo = tend.dTo; dq_crcl = tend.dq_crcl; dq_eva = tend.dq_eva; dq_rain = tend.dq_rain

        @turbo for j in 1:ydim
            for i in 1:xdim
                ts0 = Ts[i, j] + dT_ocean[i, j] + Δt * @surface_flux(i, j) / cap_surf[i, j]
                tfc = (Tc[i, j] - ts0) * cap_surf[i, j] / Δt
                ts0 = ts0 + tfc * Δt / cap_surf[i, j]
                TFc[i, j] = tfc
                Ts0_buf[i, j] = ts0

                Ta0_buf[i, j] = Ta[i, j] + dTa_crcl[i, j] + ΔT_air_factor * @atmosphere_flux(i, j)

                to0 = To[i, j] + dTo[i, j]
                tofc = Toc[i, j] - to0
                to0 = to0 + tofc
                ToFc[i, j] = tofc
                To0_buf[i, j] = to0

                q0 = q[i, j] + dq_crcl[i, j] + Δt * (dq_eva[i, j] + dq_rain[i, j])
                qfc = qc[i, j] - q0
                q0 = q0 + qfc
                qFc[i, j] = qfc
                q0_buf[i, j] = q0
            end
        end

        # Sea ice (updates cap_surf in place)
        seaice!(fields, ws.Ts0, slice, r.config.processes)

        surf = SurfaceState(ws.Ts0, ws.Ta0, ws.To0, ws.q0)
        diagnostics!(t, surf, state)

        # Advance state
        @. Ts = ws.Ts0
        @. Ta = ws.Ta0
        @. To = ws.To0
        @. q = ws.q0
    end
    return nothing
end

# A fresh copy of a state's arrays
_copy_surface(s::SurfaceState) = SurfaceState(copy(s.Ts), copy(s.Ta), copy(s.To), copy(s.q))

# The corrections of a `Stored` run. With a dataset directory they are read from its file
# unless `fields` already holds them from that file; without one, the ones in `fields`
# are used (`load_climatology` loads them)
function _use_stored_corrections!(fields::ClimateFields, jld2_dir::AbstractString, in_place::Bool)
    isempty(jld2_dir) && return nothing
    file = _flux_file(jld2_dir)
    isfile(file) || throw(ArgumentError("Stored() found no flux corrections at $file"))
    if !in_place
        @info "Loading the stored flux corrections"
        load_flux_corrections!(String(jld2_dir), fields)
    end
    return nothing
end

# The control run: `years` years from a copy of `surf_ini`. Returns the monthly records and
# their times
function _control_run!(r::ResolvedConfig, fields::ClimateFields, surf_ini::SurfaceState, CO2_ctrl, years,
    state::ModelState, ws::ModelWorkspace, acc::MonthlyAccumulator;
    ws_a::ModelWorkspace, ws_q::ModelWorkspace, observer)
    reset!(acc)   # drop what the spin-up accumulated
    @info "Control run: CO2 = $CO2_ctrl ppm, $years yr"

    surf = _copy_surface(surf_ini)
    output = MonthlyRecord[]
    times = RecordTime[]
    sizehint!(output, years * months_per_year)
    sizehint!(times, years * months_per_year)

    for t in steps(control, _control_start_year, years)
        time_loop!(t, CO2_ctrl, surf.Ts, surf.Ta, surf.To, surf.q, output, times, fields, state, ws, acc, r;
            ws_a=ws_a, ws_q=ws_q, observer=observer)
    end
    return output, times
end

# The scenario run: `years` years from a copy of `surf_ini`, after the control run `ctrl_output`
# has set the ice cover the surface masks need. Returns the monthly records and their times
function _scenario_run!(r::ResolvedConfig, fields::ClimateFields, surf_ini::SurfaceState, CO2_ctrl, years,
    ctrl_output, state::ModelState, ws::ModelWorkspace, acc::MonthlyAccumulator;
    ws_a::ModelWorkspace, ws_q::ModelWorkspace, observer)
    s = r.config.scenario

    # Regional CO2 masks belong to the scenario: the spin-up and the control
    # ran on the full CO2 everywhere
    ice_forcing = ice_climatology(ctrl_output)
    apply_co2_mask!(s.co2_mask, fields)
    apply_surface_mask!(s.co2_mask, fields, ice_forcing)

    years > 0 && @info "Scenario run: $years yr"

    # Solar-table scenarios: swap in the alternate insolation
    r.solar_table === nothing || (fields.sw_solar .= r.solar_table)

    # Forced-boundary scenarios
    _add_boundary_anomaly!(s.surface, fields)

    surf = _copy_surface(surf_ini)
    reset!(acc)

    output = MonthlyRecord[]
    times = RecordTime[]
    sizehint!(output, years * months_per_year)
    sizehint!(times, years * months_per_year)

    for t in steps(scenario, s.start_year, years)
        forcing_result = forcing(t, r)
        CO2 = forcing_result.CO2
        state.sw_solar_forcing = forcing_result.sw_solar_forcing

        # Forced-boundary experiments: overwrite Ts with climatology
        if s.surface isa BoundaryAnomaly
            surf.Ts .= clim_slice(fields.Ts_clim, data_slice(t))
        end

        # Ocean surface held at climatology plus an offset, CO2 at control
        if s.surface isa SSTOffset
            CO2 = CO2_ctrl
            Ts_clim = clim_slice(fields.Ts_clim, data_slice(t))
            @. surf.Ts = ifelse(!is_land(fields.z_topo), Ts_clim + s.surface.offset, surf.Ts)
        end

        time_loop!(t, CO2, surf.Ts, surf.Ta, surf.To, surf.q, output, times, fields, state, ws, acc, r;
            ws_a=ws_a, ws_q=ws_q, observer=observer)
    end
    return output, times
end

"""
    greb_model!(run::RunSpec, config::Config; jld2_dir="", fields=ClimateFields(),
                allow_uninitialized=false, observer=nothing)
    greb_model!(run::RunSpec, r::ResolvedConfig; ...)

Run `config`: its flux corrections (a [`SpinUp`](@ref), the [`Stored`](@ref)
ones, or none), a control run of `run.ctrl` years and a scenario run of
`run.scnr` years. Returns `(ctrl, scnr, ctrl_time, scnr_time, scnr_anomaly)`:
`ctrl` and `scnr` are vectors of [`MonthlyRecord`](@ref), `ctrl_time` and
`scnr_time` hold the [`RecordTime`](@ref) (year and month) of each record.
`scnr_anomaly` is `true` when `scnr` is the anomaly against the control's final
year: the scenario's `output` is `:anomaly` and both runs have records. It is
`false` when `scnr` holds absolute values.

| Keyword | Meaning |
|:--------|:--------|
| `jld2_dir` | The dataset directory. A `Config` is [`resolve`](@ref)d from it, and the run reads stored corrections and anomaly fields from it |
| `fields` | The loaded [`ClimateFields`](@ref) (from [`load_climatology`](@ref)). The run restores the input fields it changes when it returns, so one instance can be passed to several runs |
| `allow_uninitialized` | Accept all-zero `fields`; for precompilation and tests |
| `observer` | A function `observer(point, view)`, called twice per control and scenario step |

The observer is called at `:after_tendencies` (the step's flows are known,
the state is still the old one) and at `:after_step` (the state is updated).
`view` is a `NamedTuple` with `time` (the [`ModelTime`](@ref)) and, read from
it, `phase` (`GREBClimate.control` or `GREBClimate.scenario`), `step`, `year`
and `step_of_year`; `CO2`, the state `Ts`, `Ta`, `To`, `q`, the flows `tend` (what
[`tendencies!`](@ref) returns), `fields` and `config`. These are the model's
own arrays: read or copy them, do not write to them. The spin-up does not
call the observer. `GREBClimate.BudgetCheck` is one.
"""
greb_model!(run::RunSpec, config::Config; jld2_dir::AbstractString="", kwargs...) =
    greb_model!(run, resolve(config; jld2_dir); jld2_dir, kwargs...)

function greb_model!(run::RunSpec, r::ResolvedConfig;
    jld2_dir::AbstractString="", fields::ClimateFields=ClimateFields(),
    allow_uninitialized::Bool=false, observer=nothing)
    if !fields.loaded && !allow_uninitialized
        error("""
              greb_model! was given an uninitialized ClimateFields (all-zero climatology).

              Load the data and pass it through:
                  fields = load_climatology(jld2_dir; dataset=:ncep)
                  greb_model!(run, config; jld2_dir=jld2_dir, fields=fields)

              If a data-free run is intended (precompilation, config/scenario-plumbing
              tests), opt in explicitly with `allow_uninitialized=true`.
              """)
    end
    config = r.config
    s = config.scenario
    run.ctrl == 0 && run.scnr > 0 && s.co2_mask isa SurfaceMask &&
        throw(ArgumentError("a SurfaceMask is built from the control run's ice cover: use ctrl >= 1"))
    # The run overwrites these in place; restore them so a reused `fields`
    # does not carry one run's changes into the next
    # Stored corrections that `fields` already holds from this file are neither
    # reloaded nor saved: a `Stored` run does not write them
    saved_flux_source = fields.flux_source
    flux_file = isempty(jld2_dir) ? "" : _flux_file(jld2_dir)
    flux_in_place = config.corrections isa Stored && !isempty(flux_file) && fields.flux_source == flux_file
    saved = [name => copy(getfield(fields, name)) for name in _mutated_fields(config)
             if !(flux_in_place && name in (:Ts_flux_correction, :q_flux_correction, :To_flux_correction))]
    try
        state = ModelState()

        # ── 1. Initialisation ───────────────────────────────────
        # The anomaly files are read once per `fields`, directory and source
        if s.surface isa BoundaryAnomaly
            dir, source = String(jld2_dir), s.surface.source
            fields.boundary_anomaly_source == (dir, source) || load_boundary_anomaly!(dir, fields, source)
        end

        ini = init_model!(r, fields)
        surf_ini, CO2_ctrl = ini.surf, ini.CO2_ctrl

        # Workspace and accumulator. `ws_a`/`ws_q` are separate `circulation!`
        # scratch spaces so the Ta/q circulation calls inside `tendencies!` can
        # run concurrently on `Threads.nthreads() > 1`
        ws = ModelWorkspace()
        ws_a = ModelWorkspace()
        ws_q = ModelWorkspace()
        acc = MonthlyAccumulator()

        # ── 2. Flux-correction spin-up ──────────────────────────
        corrections = config.corrections
        if corrections isa Stored
            _use_stored_corrections!(fields, jld2_dir, flux_in_place)
        elseif corrections isa SpinUp
            @info "Flux-correction spin-up: CO2 = $CO2_ctrl ppm, $(corrections.years) yr"
            qflux_correction!(CO2_ctrl, surf_ini.Ts, surf_ini.Ta, surf_ini.To, surf_ini.q, fields, state, r, ws,
                corrections.years; ws_a=ws_a, ws_q=ws_q)
        else
            @info "No flux corrections"
        end

        # ── 3. Control run and 4. scenario run ──────────────────
        ctrl_output, ctrl_time = _control_run!(r, fields, surf_ini, CO2_ctrl, run.ctrl, state, ws, acc;
            ws_a=ws_a, ws_q=ws_q, observer=observer)
        scnr_output, scnr_time = _scenario_run!(r, fields, surf_ini, CO2_ctrl, run.scnr, ctrl_output, state, ws, acc;
            ws_a=ws_a, ws_q=ws_q, observer=observer)

        # Post-processing: anomalies against the control's final year
        scnr_anomaly = s.output === :anomaly && !isempty(ctrl_output) && !isempty(scnr_output)
        if scnr_anomaly
            ctrl_clim = monthly_climatology(ctrl_output)
            scnr_output = scenario_anomalies(scnr_output, scnr_time, ctrl_clim)
        end

        return (ctrl=ctrl_output, scnr=scnr_output, ctrl_time=ctrl_time, scnr_time=scnr_time,
            scnr_anomaly=scnr_anomaly)
    finally
        for (name, a) in saved
            getfield(fields, name) .= a
        end
        # The arrays are back to what they held before the run, so the marker is too
        fields.flux_source = saved_flux_source
    end
end
