"""
    init_model!(r::ResolvedConfig, fields::ClimateFields)

One-time per-run setup: resets the CO2 mask to full CO2, applies the climatologies the
[`Processes`](@ref) options replace (flat topography, cloud cover, humidity,
mixed layer) and zeroes the flux corrections for `NoCorrections()`, derives
the fields that follow from them ([`derive_fields!`](@ref)), then computes the
control-run initial state. Returns
`(Ts_ini, Ta_ini, To_ini, q_ini, CO2_ctrl)`.
"""
function init_model!(r::ResolvedConfig, fields::ClimateFields)
    p = r.config.processes

    # The spin-up and the control run on the full CO2 everywhere; a regional
    # mask is set at the start of the scenario. The reset also clears what an
    # earlier regional run left in a reused `fields`
    fields.co2_part .= 1.0f0

    Ts_clim = fields.Ts_clim
    z_topo = fields.z_topo

    # ── Sensitivity experiment overrides ─────────────────────────
    if p.topography === :flat
        @. z_topo = min(z_topo, 1.0f0)       # constant topography
    end

    if p.clouds === :none
        fields.cloud_clim .= 0.0f0  # zero cloud climatology
    end

    if r.config.corrections isa NoCorrections
        fields.Ts_flux_correction .= 0.0f0
        fields.q_flux_correction .= 0.0f0
        fields.To_flux_correction .= 0.0f0
    end

    # Climatology modifications
    if p.hydrology === :none
        fields.q_clim .= 0.0f0  # zero out humidity climatology
    end

    if p.clouds === :uniform
        fields.cloud_clim .= 0.7f0           # constant cloud cover (2xCO2 deconstruction)
    end

    if p.humidity === :uniform
        fields.q_clim .= 0.0052f0          # constant water vapor
    end

    if p.ocean === :mixed_layer
        fields.mld_clim .= d_ocean       # no deep ocean
    end

    derive_fields!(fields, p)

    # ── Initial conditions from last time step of climatology ────
    Ts_ini = Ts_clim[:, :, nstep_yr] |> copy          # surface temperature
    Ta_ini = copy(Ts_ini)                           # air temperature = Tsurf
    To_ini = fields.To_clim[:, :, nstep_yr] |> copy  # deep ocean temperature
    q_ini = fields.q_clim[:, :, nstep_yr] |> copy    # atmospheric water vapor

    # ── Control CO2 level ───────────────────────────────────────
    CO2_ctrl = p.co2 ? r.config.scenario.control_co2 : 0.0f0

    return (Ts_ini=Ts_ini, Ta_ini=Ta_ini, To_ini=To_ini,
        q_ini=q_ini, CO2_ctrl=CO2_ctrl)
end

"The climatologies a [`BoundaryAnomaly`](@ref) adds its anomalies to at scenario start."
const _BOUNDARY_FIELDS = (:Ts_clim, :u_clim, :v_clim, :omega_clim, :wind_speed_clim)

_add_boundary_anomaly!(::SurfaceForcing, fields::ClimateFields) = fields

function _add_boundary_anomaly!(b::BoundaryAnomaly, fields::ClimateFields)
    if b.source === :cmip5_rcp85
        @info "Applying CMIP5 RCP8.5 climate change forcing"
        suffix = :_anom_cc
    else
        @info "Applying ERA-Interim $(b.source === :elnino ? "El Niño" : "La Niña") forcing"
        suffix = :_anom_enso
    end
    for name in _BOUNDARY_FIELDS
        getfield(fields, name) .+= getfield(fields, Symbol(name, suffix))
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
    p = c.processes
    names = [:sw_solar, :Ts_flux_correction, :q_flux_correction, :To_flux_correction]
    c.scenario.surface isa BoundaryAnomaly && append!(names, _BOUNDARY_FIELDS)
    p.topography === :flat && push!(names, :z_topo)
    p.clouds === :observed || push!(names, :cloud_clim)
    (p.hydrology === :none || p.humidity === :uniform) && push!(names, :q_clim)
    p.ocean === :mixed_layer && push!(names, :mld_clim)
    return names
end

# The control run and the spin-up carry this year as their label
const _control_start_year = 1970

"""
    qflux_correction!(CO2_ctrl, Ts, Ta, q, To, fields, state, r::ResolvedConfig, ws, years; ws_a=ws, ws_q=ws)

Runs `years` years of `tendencies!` to derive the ocean/atmosphere flux
corrections (`fields.Ts_flux_correction`/`q_flux_correction`/`To_flux_correction`) that make the
control climate match observed climatology. Mutates `Ts`/`Ta`/`q`/`To` in
place as it integrates. `ws_a`/`ws_q` are forwarded to [`tendencies!`](@ref).
"""
function qflux_correction!(CO2_ctrl, Ts, Ta, q, To, fields::ClimateFields, state::ModelState, r::ResolvedConfig, ws::ModelWorkspace, years;
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
        seaice!(ws.Ts0, fields, slice, r.config.processes)

        surf = SurfaceState(ws.Ts0, ws.Ta0, ws.To0, ws.q0)
        diagnostics!(t, surf, state)

        # Advance state
        @. Ts = ws.Ts0
        @. Ta = ws.Ta0
        @. q = ws.q0
        @. To = ws.To0
    end
    return nothing
end

"""
    greb_model!(run::RunSpec, config::Config; jld2_dir="", fields=ClimateFields(),
                allow_uninitialized=false, observer=nothing)
    greb_model!(run::RunSpec, r::ResolvedConfig; ...)

Run `config`: its flux corrections (a [`SpinUp`](@ref), the [`Stored`](@ref)
ones, or none), a control run of `run.ctrl` years and a scenario run of
`run.scnr` years. Returns `(ctrl, scnr, ctrl_time, scnr_time)`: `ctrl` and
`scnr` are vectors of [`MonthlyRecord`](@ref), `ctrl_time` and `scnr_time` hold
the [`RecordTime`](@ref) (year and month) of each record. `scnr` is the
anomaly against the control's final year when the scenario's `output` is
`:anomaly`.

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
    time_ctrl, time_scnr = run.ctrl, run.scnr
    time_ctrl == 0 && time_scnr > 0 && s.co2_mask isa SurfaceMask &&
        throw(ArgumentError("a SurfaceMask is built from the control run's ice cover: use ctrl >= 1"))
    is_forced_boundary = s.surface isa BoundaryAnomaly
    # The run overwrites these in place; restore them so a reused `fields`
    # does not carry one run's changes into the next
    saved = [name => copy(getfield(fields, name)) for name in _mutated_fields(config)]
    try

    state = ModelState()

    # ── 1. Initialisation ───────────────────────────────────────
    # The anomaly files are read once per `fields`, directory and source
    if is_forced_boundary
        dir, source = String(jld2_dir), s.surface.source
        loaded = source === :cmip5_rcp85 ? fields.anom_cc_source == dir :
                 fields.anom_enso_source == (dir, source)
        loaded || load_boundary_anomaly!(dir, fields, source)
    end

    ini = init_model!(r, fields)
    Ts_ini = ini.Ts_ini;
    Ta_ini = ini.Ta_ini
    To_ini = ini.To_ini;
    q_ini = ini.q_ini
    CO2_ctrl = ini.CO2_ctrl

    # Workspace and accumulator. `ws_a`/`ws_q` are separate `circulation!`
    # scratch spaces so the Ta/q circulation calls inside `tendencies!` can
    # run concurrently on `Threads.nthreads() > 1`
    ws = ModelWorkspace()
    ws_a = ModelWorkspace()
    ws_q = ModelWorkspace()
    acc = MonthlyAccumulator()

    # ── 2. Flux-correction spin-up ──────────────────────────────
    corrections = config.corrections
    if corrections isa Stored
        # Without a directory the corrections already in `fields` are used
        # (`load_climatology` loads them)
        if !isempty(jld2_dir)
            file = joinpath(jld2_dir, "climatology", "flux_corrections.jld2")
            isfile(file) || throw(ArgumentError("Stored() found no flux corrections at $file"))
            @info "Loading the stored flux corrections"
            load_flux_corrections!(String(jld2_dir), fields)
        end
    elseif corrections isa SpinUp
        @info "Flux-correction spin-up: CO2 = $CO2_ctrl ppm, $(corrections.years) yr"
        qflux_correction!(CO2_ctrl, Ts_ini, Ta_ini, q_ini, To_ini, fields, state, r, ws, corrections.years;
            ws_a=ws_a, ws_q=ws_q)
    else
        @info "No flux corrections"
    end

    # Reset accumulators after spin-up
    reset!(acc)

    # ── 3. Control run ──────────────────────────────────────────
    @info "Control run: CO2 = $CO2_ctrl ppm, $time_ctrl yr"

    # Initialize state arrays
    Ts = copy(Ts_ini);
    Ta = copy(Ta_ini)
    To = copy(To_ini);
    q = copy(q_ini)
    state.sw_solar_forcing = 1.0f0

    ctrl_output = MonthlyRecord[]
    ctrl_time = RecordTime[]
    sizehint!(ctrl_output, time_ctrl * months_per_year)
    sizehint!(ctrl_time, time_ctrl * months_per_year)

    for t in steps(control, _control_start_year, time_ctrl)
        time_loop!(t, CO2_ctrl, Ts, Ta, q, To, ctrl_output, ctrl_time, fields, state, ws, acc, r;
            ws_a=ws_a, ws_q=ws_q, observer=observer)
    end

    # ── Build ice climatology from control output ───────────────
    ice_forcing = ice_climatology(ctrl_output)

    # Regional CO2 masks belong to the scenario: the spin-up and the control
    # above ran on the full CO2 everywhere
    apply_co2_mask!(s.co2_mask, fields)
    apply_surface_mask!(s.co2_mask, fields, ice_forcing)

    # ── 4. Scenario run ─────────────────────────────────────────
    time_scnr > 0 && @info "Scenario run: $time_scnr yr"

    # Solar-table scenarios: swap in the alternate insolation
    r.solar_table === nothing || (fields.sw_solar .= r.solar_table)

    # Forced-boundary scenarios
    _add_boundary_anomaly!(s.surface, fields)

    # Reset state to initial conditions
    Ts .= Ts_ini;
    Ta .= Ta_ini
    q .= q_ini;
    To .= To_ini

    state.sw_solar_forcing = 1.0f0
    reset!(acc)

    scnr_output = MonthlyRecord[]
    scnr_time = RecordTime[]
    sizehint!(scnr_output, time_scnr * months_per_year)
    sizehint!(scnr_time, time_scnr * months_per_year)

    for t in steps(scenario, s.start_year, time_scnr)
        forcing_result = forcing(t, r)
        CO2 = forcing_result.CO2
        state.sw_solar_forcing = forcing_result.sw_solar_forcing

        # Forced‑boundary experiments: overwrite Ts with climatology
        if is_forced_boundary
            Ts .= clim_slice(fields.Ts_clim, data_slice(t))
        end

        # Ocean surface held at climatology plus an offset, CO2 at control
        if s.surface isa SSTOffset
            CO2 = CO2_ctrl
            Ts_clim = clim_slice(fields.Ts_clim, data_slice(t))
            @. Ts = ifelse(!is_land(fields.z_topo), Ts_clim + s.surface.offset, Ts)
        end

        time_loop!(t, CO2, Ts, Ta, q, To, scnr_output, scnr_time, fields, state, ws, acc, r;
            ws_a=ws_a, ws_q=ws_q, observer=observer)
    end

    # Post‑processing: anomalies against the control's final year
    if s.output === :anomaly && !isempty(ctrl_output) && !isempty(scnr_output)
        ctrl_clim = monthly_climatology(ctrl_output)
        scnr_output = scenario_anomalies(scnr_output, scnr_time, ctrl_clim)
    end

    return (ctrl=ctrl_output, scnr=scnr_output, ctrl_time=ctrl_time, scnr_time=scnr_time)

    finally
        for (name, a) in saved
            getfield(fields, name) .= a
        end
    end
end
