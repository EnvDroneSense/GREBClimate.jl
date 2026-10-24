"""
    tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, slice, r::ResolvedConfig; ws_a=ws, ws_q=ws)

Run one timestep's physics pipeline and return a named tuple of every
intermediate flux/tendency needed by [`diagnostics!`](@ref) and the caller's own
state update. [`circulation!`](@ref) (temperature, then humidity) starts first,
then [`SWradiation!`](@ref) → [`LWradiation!`](@ref) → sensible heat →
[`hydro!`](@ref) → [`deep_ocean!`](@ref) run on `ws`.

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
        t_a = Threads.@spawn circulation!(ws_a.dTa_crcl, Ta, z_air, fields, ws_a, slice, p)
        t_q = Threads.@spawn circulation!(ws_q.dq_crcl, q, z_vapor, fields, ws_q, slice, p)
    else
        circulation!(ws_a.dTa_crcl, Ta, z_air, fields, ws_a, slice, p)
        circulation!(ws_q.dq_crcl, q, z_vapor, fields, ws_q, slice, p)
    end

    # Short-wave radiation -> albedo, SW flux
    sw_out = SWradiation!(ws, Ts, fields, state, slice, p)

    # Long-wave radiation -> LW_surf, LW_up, LW_down, emissivity
    lw_out = LWradiation!(ws, Ts, Ta, q, CO2, fields, slice, p)

    # Sensible heat flux
    Q_sens = ws.Q_sens
    if p.atmosphere
        @. Q_sens = ct_sens * (Ta - Ts)
    else
        fill!(Q_sens, 0.0f0)
    end

    # Hydrological cycle -> latent heat + evaporation/rain tendencies
    hy_out = hydro!(ws, Ts, q, fields, slice, p, r.hydrology)

    # Deep ocean coupling
    do_out = deep_ocean!(ws, Ts, To, fields, slice, p)

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
