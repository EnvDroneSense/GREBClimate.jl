"""
    tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, slice, r::ResolvedConfig; ws_a=ws, ws_q=ws)

Runs one timestep's physics pipeline: [`SWradiation!`](@ref) →
[`LWradiation!`](@ref) → sensible heat → [`hydro!`](@ref) →
[`circulation!`](@ref) (temperature, then humidity) → [`deep_ocean!`](@ref)
and returns a named tuple of every intermediate flux/tendency needed by
[`diagnostics!`](@ref) and the caller's own state update.

The two `circulation!` calls are independent of each other and of every
other stage (each reads only pre-timestep state and writes disjoint
buffers), so when the caller supplies distinct `ws_a`/`ws_q` workspaces
*and* `Threads.nthreads() > 1`, they run concurrently via `Threads.@spawn`
while the remaining stages run on `ws`. With the default `ws_a=ws_q=ws` they
run one after the other.
"""
function tendencies!(CO2, Ts, Ta, To, q, fields::ClimateFields, state::ModelState, ws::ModelWorkspace,
    slice::Int, r::ResolvedConfig; ws_a::ModelWorkspace=ws, ws_q::ModelWorkspace=ws)

    p = r.config.processes
    parallel = Threads.nthreads() > 1 && ws_a !== ws_q

    # Atmospheric circulation - temperature/water-vapor diffusion/advection.
    if parallel
        t_a = Threads.@spawn circulation!(Ta, z_air, ws_a.dTa_crcl, fields, ws_a, slice, p)
        t_q = Threads.@spawn circulation!(q, z_vapor, ws_q.dq_crcl, fields, ws_q, slice, p)
    else
        circulation!(Ta, z_air, ws_a.dTa_crcl, fields, ws_a, slice, p)
        circulation!(q, z_vapor, ws_q.dq_crcl, fields, ws_q, slice, p)
    end

    # Short-wave radiation -> albedo, SW flux
    sw_out = SWradiation!(Ts, fields, state, slice, p, ws)

    # Long-wave radiation -> LW_surf, LW_up, LW_down, emissivity
    lw_out = LWradiation!(Ts, Ta, q, CO2, fields, slice, p, ws)

    # Sensible heat flux
    Q_sens = ws.Q_sens
    if p.atmosphere
        @. Q_sens = ct_sens * (Ta - Ts)
    else
        fill!(Q_sens, 0.0f0)
    end

    # Hydrological cycle -> latent heat + evaporation/rain tendencies
    hy_out = hydro!(Ts, q, fields, slice, p, r.hydrology, ws)

    # Deep ocean coupling
    do_out = deep_ocean!(Ts, To, fields, slice, p, ws)

    if parallel
        wait(t_a)
        wait(t_q)
    end

    return (albedo=sw_out.albedo,
        SW=sw_out.SW,
        ice_cover=sw_out.ice_cover,
        LW_surf=lw_out.LW_surf,
        Q_lat=hy_out.Q_lat,
        Q_sens=Q_sens,
        Q_lat_air=hy_out.Q_lat_air,
        dq_eva=hy_out.dq_eva,
        dq_rain=hy_out.dq_rain,
        dq_crcl=ws_q.dq_crcl,
        dTa_crcl=ws_a.dTa_crcl,
        dT_ocean=do_out.dT_ocean,
        dTo=do_out.dTo,
        LW_down=lw_out.LW_down,
        LW_up=lw_out.LW_up,
        em=lw_out.em)
end
