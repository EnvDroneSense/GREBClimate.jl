# ── Ghost-cell helpers ───────────────────────────────────────────────
# Ghosted buffer layout and why it exists: see `nghost` in `constants.jl`.

"Refresh the wrap-around ghost rows of a `(xghost, ydim)` buffer."
@inline function refresh_ghosts!(P::AbstractMatrix{Float32})
    @inbounds for k in axes(P, 2), h in 1:nghost
        P[h, k] = P[xdim+h, k]              # west ghosts <- east edge
        P[xdim+nghost+h, k] = P[nghost+h, k]    # east ghosts <- west edge
    end
    return nothing
end

"Refresh the wrap-around ghost entries of a length-`xghost` vector."
@inline function refresh_ghosts!(v::AbstractVector{Float32})
    @inbounds for h in 1:nghost
        v[h] = v[xdim+h]
        v[xdim+nghost+h] = v[nghost+h]
    end
    return nothing
end

"Copy an `(xdim, ydim)` field into the ghosted buffer `P` and fill its ghosts."
function to_ghosted!(P::Matrix{Float32}, A::Matrix{Float32})
    @inbounds for k in 1:ydim
        copyto!(P, (k - 1) * xghost + nghost + 1, A, (k - 1) * xdim + 1, xdim)
    end
    refresh_ghosts!(P)
    return P
end

function to_ghosted!(P::AbstractMatrix{Float32}, A::AbstractMatrix{<:Real})
    @inbounds for k in 1:ydim
        @simd for i in 1:xdim
            P[i+nghost, k] = A[i, k]
        end
    end
    refresh_ghosts!(P)
    return P
end

"""
    convergence!(T1, fields::ClimateFields, timestate, ws::CirculationWorkspace)

Moisture flux convergence from `T1` (specific humidity, `[kg/kg]`) and the
current `fields.omegaclim` (vertical velocity), writing the tendency into
`ws.dX_conv`. Implements Eq. 18 from Stassen et al. (2019).
"""
function convergence!(T1, fields::ClimateFields, timestate, ws::CirculationWorkspace)
    omega = @view fields.omegaclim[:, :, timestate.ityr]

    @. ws.dX_conv = -T1 * omega * const_factor
    return nothing
end

# Same tendency, reading a ghosted field.
function _convergence!(Xp::Matrix{Float32}, fields::ClimateFields, timestate, ws::CirculationWorkspace)
    omega = @view fields.omegaclim[:, :, timestate.ityr]
    dX_conv = ws.dX_conv

    @inbounds for k in 1:ydim
        @turbo for i in 1:xdim
            dX_conv[i, k] = -Xp[i+nghost, k] * omega[i, k] * const_factor
        end
    end
    return nothing
end

"Select `wz_air`/`wz_vapor` for a scale height, erroring on anything else."
@inline function _wz_for(h_scl, fields::ClimateFields)
    if h_scl == z_air
        return fields.wz_air
    elseif h_scl == z_vapor
        return fields.wz_vapor
    else
        error("Invalid h_scl = $h_scl (must be z_air or z_vapor)")
    end
end

"""
    diffusion!(T1, h_scl, fields::ClimateFields, ws::CirculationWorkspace, timestate)

Meridional + zonal diffusion of `T1` (temperature or humidity), writing the
tendency into `ws.dX_diff`. `h_scl` (`z_air` or `z_vapor`) selects the
topographic weighting field.
"""
function diffusion!(T1, h_scl, fields::ClimateFields, ws::CirculationWorkspace, timestate)
    wz = _wz_for(h_scl, fields)
    to_ghosted!(ws.X_work, T1)
    to_ghosted!(ws.wz_ghost, wz)
    _diffusion!(ws.X_work, ws.wz_ghost, ws)
    return nothing
end

# Core kernel. `Tp`/`wzp` are ghosted with valid ghost rows; the tendency lands in
# the plain `(xdim, ydim)` `ws.dX_diff`.
function _diffusion!(Tp::Matrix{Float32}, wzp::Matrix{Float32}, ws::CirculationWorkspace)
    dX_diff = ws.dX_diff

    # Precomputed geometry/coefficients
    ccy = ccy_diff
    ccx = ccx_diff
    is_polar = IS_POLAR

    term_north = ws.term_north
    term_south = ws.term_south
    T1h = ws.T1h
    dTxh = ws.dTxh

    # ----- Precompute k-independent terms for the poles -----
    @turbo for i in 1:xdim
        ip = i + nghost
        # For k == 1 (North Pole)
        term_north[i] = ccy * wzp[ip, 2] * (Tp[ip, 2] - Tp[ip, 1])
        # For k == ydim (South Pole)
        term_south[i] = ccy * wzp[ip, ydim-1] * (Tp[ip, ydim-1] - Tp[ip, ydim])
    end

    @inbounds for k in 1:ydim
        # ----- Meridional diffusion -----
        if k == 1
            @turbo for i in 1:xdim
                dX_diff[i, k] = wzp[i+nghost, k] * term_north[i]
            end
        elseif k == ydim
            @turbo for i in 1:xdim
                dX_diff[i, k] = wzp[i+nghost, k] * term_south[i]
            end
        else
            # Mid-latitudes: no precomputation possible (depends on k-1, k+1)
            @turbo for i in 1:xdim
                dX_diff[i, k] = wzp[i+nghost, k] * ccy * (
                    wzp[i+nghost, k-1] * (Tp[i+nghost, k-1] - Tp[i+nghost, k]) +
                    wzp[i+nghost, k+1] * (Tp[i+nghost, k+1] - Tp[i+nghost, k])
                )
            end
        end

        # ----- Zonal diffusion -----
        if !is_polar[k]   # mid-latitudes, normal time step
            cc = ccx[k] * 0.05f0
            @turbo for j in 1:xdim
                # Ghost offsets: j-3 -> j, j-2 -> j+1, j-1 -> j+2, j -> j+3, ...
                wm3 = wzp[j, k];   wm2 = wzp[j+1, k]; wm1 = wzp[j+2, k]
                w0 = wzp[j+3, k]
                wp1 = wzp[j+4, k]; wp2 = wzp[j+5, k]; wp3 = wzp[j+6, k]
                tm3 = Tp[j, k];    tm2 = Tp[j+1, k];  tm1 = Tp[j+2, k]
                t0 = Tp[j+3, k]
                tp1 = Tp[j+4, k];  tp2 = Tp[j+5, k];  tp3 = Tp[j+6, k]

                dTx = cc * (
                    10.0f0 * (wm1 * (tm1 - t0) + wp1 * (tp1 - t0)) +
                    4.0f0 * (wm2 * (tm2 - tm1) + wm1 * (t0 - tm1)) +
                    4.0f0 * (wp1 * (t0 - tp1) + wp2 * (tp2 - tp1)) +
                    1.0f0 * (wm3 * (tm3 - tm2) + wm2 * (tm1 - tm2)) +
                    1.0f0 * (wp2 * (tp1 - tp2) + wp3 * (tp3 - tp2))
                )
                dX_diff[j, k] += w0 * dTx
            end
        else   # polar regions - sub-timestepping
            # Number of sub-steps for stability (precomputed, depends only on k)
            time2 = POLAR_DIFF_TIME2[k]
            cc2 = POLAR_DIFF_CCX2[k] * 0.05f0

            # Copy current row (ghosts included) into the temporary buffer
            @simd for i in 1:xghost
                T1h[i] = Tp[i, k]
            end

            for _ in 1:time2
                # Jacobi
                @turbo for j in 1:xdim
                    wm3 = wzp[j, k];   wm2 = wzp[j+1, k]; wm1 = wzp[j+2, k]
                    wp1 = wzp[j+4, k]; wp2 = wzp[j+5, k]; wp3 = wzp[j+6, k]
                    tm3 = T1h[j];      tm2 = T1h[j+1];    tm1 = T1h[j+2]
                    t0 = T1h[j+3]
                    tp1 = T1h[j+4];    tp2 = T1h[j+5];    tp3 = T1h[j+6]

                    dTxh[j] = cc2 * (
                        10.0f0 * (wm1 * (tm1 - t0) + wp1 * (tp1 - t0)) +
                        4.0f0 * (wm2 * (tm2 - tm1) + wm1 * (t0 - tm1)) +
                        4.0f0 * (wp1 * (t0 - tp1) + wp2 * (tp2 - tp1)) +
                        1.0f0 * (wm3 * (tm3 - tm2) + wm2 * (tm1 - tm2)) +
                        1.0f0 * (wp2 * (tp1 - tp2) + wp3 * (tp3 - tp2))
                    )
                end
                @turbo for j in 1:xdim
                    t0 = T1h[j+3]
                    dq = ifelse(dTxh[j] <= -t0, -0.9f0 * t0, dTxh[j])
                    T1h[j+3] = t0 + dq
                end
                refresh_ghosts!(T1h)
            end

            # Add total change (scaled by outer wz) to output buffer
            @turbo for i in 1:xdim
                dX_diff[i, k] += wzp[i+nghost, k] * (T1h[i+nghost] - Tp[i+nghost, k])
            end
        end
    end

    return nothing
end

"""
    advection!(T1, h_scl, fields::ClimateFields, ws::CirculationWorkspace, timestate, cfg::PhysicsConfig)

Meridional + zonal advection of `T1` (temperature or humidity), writing the
tendency into `ws.dX_adv`. Gated by `cfg.log_hadv`/`cfg.log_vadv` depending on
`h_scl`.
"""
function advection!(T1, h_scl, fields::ClimateFields, ws::CirculationWorkspace, timestate, cfg::PhysicsConfig)
    # Disable advection for water vapour or heat according to switches
    if (h_scl == z_vapor && !cfg.log_vadv) || (h_scl == z_air && !cfg.log_hadv)
        fill!(ws.dX_adv, 0.0f0)
        return nothing
    end
    wz = _wz_for(h_scl, fields)
    to_ghosted!(ws.X_work, T1)
    to_ghosted!(ws.wz_ghost, wz)
    _advection!(ws.X_work, ws.wz_ghost, fields, ws, timestate)
    return nothing
end

# Core kernel. `Tp`/`wzp` are ghosted; the switch check has already run.
function _advection!(Tp::Matrix{Float32}, wzp::Matrix{Float32}, fields::ClimateFields,
                     ws::CirculationWorkspace, timestate)
    dX_adv = ws.dX_adv

    # Extract 2D views for current time step
    vclim_p_t = @view fields.vclim_p[:, :, timestate.ityr]
    vclim_m_t = @view fields.vclim_m[:, :, timestate.ityr]
    uclim_p_t = @view fields.uclim_p[:, :, timestate.ityr]
    uclim_m_t = @view fields.uclim_m[:, :, timestate.ityr]

    # Precomputed constants
    ccy = ccy_adv
    ccx = ccx_adv
    is_polar = IS_POLAR

    T1h = ws.T1h
    dTxh = ws.dTxh

    @inbounds for k in 1:ydim
        # ----- Meridional (v) advection -----
        if k == 1          # North Pole
            @turbo for j in 1:xdim
                v_p = vclim_p_t[j, k]
                dX_adv[j, k] = ccy * v_p * (
                    wzp[j+nghost, 2] * (Tp[j+nghost, 1] - Tp[j+nghost, 2]) +
                    wzp[j+nghost, 3] * (Tp[j+nghost, 1] - Tp[j+nghost, 3])
                ) / 3.0f0
            end
        elseif k == 2
            @turbo for j in 1:xdim
                v_m = vclim_m_t[j, k]
                v_p = vclim_p_t[j, k]
                dX_adv[j, k] = ccy * (
                    -v_m * wzp[j+nghost, 1] * (Tp[j+nghost, 2] - Tp[j+nghost, 1]) +
                    v_p * (wzp[j+nghost, 3] * (Tp[j+nghost, 2] - Tp[j+nghost, 3]) +
                           wzp[j+nghost, 4] * (Tp[j+nghost, 2] - Tp[j+nghost, 4])) / 3.0f0
                )
            end
        elseif k >= 3 && k <= ydim-2
            km1, km2 = k-1, k-2
            kp1, kp2 = k+1, k+2
            @turbo for j in 1:xdim
                v_m = vclim_m_t[j, k]
                v_p = vclim_p_t[j, k]
                dX_adv[j, k] = ccy * (
                    -v_m * (wzp[j+nghost, km1] * (Tp[j+nghost, k] - Tp[j+nghost, km1]) +
                            wzp[j+nghost, km2] * (Tp[j+nghost, k] - Tp[j+nghost, km2])) +
                    v_p * (wzp[j+nghost, kp1] * (Tp[j+nghost, k] - Tp[j+nghost, kp1]) +
                           wzp[j+nghost, kp2] * (Tp[j+nghost, k] - Tp[j+nghost, kp2]))
                ) / 3.0f0
            end
        elseif k == ydim-1
            km1, km2 = k-1, k-2
            kp1 = k+1
            @turbo for j in 1:xdim
                v_m = vclim_m_t[j, k]
                v_p = vclim_p_t[j, k]
                dX_adv[j, k] = ccy * (
                    -v_m * (wzp[j+nghost, km1] * (Tp[j+nghost, k] - Tp[j+nghost, km1]) +
                            wzp[j+nghost, km2] * (Tp[j+nghost, k] - Tp[j+nghost, km2])) / 3.0f0 +
                    v_p * wzp[j+nghost, kp1] * (Tp[j+nghost, k] - Tp[j+nghost, kp1])
                )
            end
        else               # k == ydim (South Pole)
            km1, km2 = k-1, k-2
            @turbo for j in 1:xdim
                v_m = vclim_m_t[j, k]
                dX_adv[j, k] = ccy * (
                    -v_m * (wzp[j+nghost, km1] * (Tp[j+nghost, k] - Tp[j+nghost, km1]) +
                            wzp[j+nghost, km2] * (Tp[j+nghost, k] - Tp[j+nghost, km2]))
                ) / 3.0f0
            end
        end

        # ----- Zonal (u) advection -----
        if !is_polar[k]   # mid-latitudes, normal timestep
            cc = ccx[k]
            @turbo for j in 1:xdim
                wm2 = wzp[j+1, k]; wm1 = wzp[j+2, k]
                wp1 = wzp[j+4, k]; wp2 = wzp[j+5, k]
                tm2 = Tp[j+1, k];  tm1 = Tp[j+2, k]
                t0 = Tp[j+3, k]
                tp1 = Tp[j+4, k];  tp2 = Tp[j+5, k]
                u_m = uclim_m_t[j, k]
                u_p = uclim_p_t[j, k]
                dX_adv[j, k] += cc * (
                    -u_m * (wm1 * (t0 - tm1) + wm2 * (t0 - tm2)) +
                    u_p * (wp1 * (t0 - tp1) + wp2 * (t0 - tp2))
                ) / 3.0f0
            end
        else # polar regions - sub-timestepping
            # Number of sub-steps (CFL stability. Precomputed, depends only on k)
            time2 = POLAR_ADV_TIME2[k]
            ccx2 = POLAR_ADV_CCX2[k]

            # Copy current row (ghosts included) into the temporary buffer
            @simd for i in 1:xghost
                T1h[i] = Tp[i, k]
            end

            for _ in 1:time2
                # Jacobi
                @turbo for j in 1:xdim
                    wm3 = wzp[j, k];   wm2 = wzp[j+1, k]; wm1 = wzp[j+2, k]
                    wp1 = wzp[j+4, k]; wp2 = wzp[j+5, k]; wp3 = wzp[j+6, k]
                    tm3 = T1h[j];      tm2 = T1h[j+1];    tm1 = T1h[j+2]
                    t0 = T1h[j+3]
                    tp1 = T1h[j+4];    tp2 = T1h[j+5];    tp3 = T1h[j+6]
                    u_m = uclim_m_t[j, k]
                    u_p = uclim_p_t[j, k]

                    dTxh[j] = ccx2 * (
                        -u_m * (10.0f0 * wm1 * (t0 - tm1) +
                                4.0f0 * wm2 * (tm1 - tm2) +
                                1.0f0 * wm3 * (tm2 - tm3)) +
                        u_p * (10.0f0 * wp1 * (t0 - tp1) +
                               4.0f0 * wp2 * (tp1 - tp2) +
                               1.0f0 * wp3 * (tp2 - tp3))
                    ) / 20.0f0
                end
                @turbo for j in 1:xdim
                    # Stability clamp (avoid negative water vapour)
                    t0 = T1h[j+3]
                    dq = ifelse(dTxh[j] <= -t0, -0.9f0 * t0, dTxh[j])
                    T1h[j+3] = t0 + dq
                end
                refresh_ghosts!(T1h)
            end

            # Add total change to the output buffer
            @turbo for i in 1:xdim
                dX_adv[i, k] += T1h[i+nghost] - Tp[i+nghost, k]
            end
        end
    end

    return nothing
end

"""
    circulation!(X_in, h_scl, dX_out, fields::ClimateFields, ws::CirculationWorkspace, timestate, cfg::PhysicsConfig)

Sub-steps `X_in` through `ntime` iterations of [`diffusion!`](@ref),
[`advection!`](@ref), and [`convergence!`](@ref) (each gated by the relevant
`cfg.log_*` switch), writing the total change into `dX_out`. The sub-step
loop is a genuine sequential recurrence and is not parallelized.
"""
function circulation!(X_in, h_scl, dX_out, fields::ClimateFields, ws::CirculationWorkspace, timestate, cfg::PhysicsConfig)
    # Early exit if atmospheric processes disabled
    if (!cfg.log_atmos_dmc || !cfg.log_crcl_dmc || !cfg.log_crcl_drsp)
        fill!(dX_out, 0.0f0)
        return nothing
    end

    # Precompute flags
    do_diff_v = cfg.log_vdif && h_scl == z_vapor
    do_diff_h = cfg.log_hdif && h_scl == z_air
    do_adv_v = cfg.log_vadv && h_scl == z_vapor
    do_adv_h = cfg.log_hadv && h_scl == z_air
    do_conv = cfg.log_conv && h_scl == z_vapor

    # `wz` is static for the whole run, so its ghosted copy is built once per
    # call and reused across all `ntime` sub-steps.
    wzp = ws.wz_ghost
    to_ghosted!(wzp, _wz_for(h_scl, fields))

    Xp = ws.X_work
    to_ghosted!(Xp, X_in)

    dX_diff = ws.dX_diff
    dX_adv = ws.dX_adv
    dX_conv = ws.dX_conv
    fill!(dX_diff, 0.0f0)
    fill!(dX_adv, 0.0f0)
    fill!(dX_conv, 0.0f0)

    for _tt in 1:ntime
        (do_diff_v || do_diff_h) && _diffusion!(Xp, wzp, ws)
        (do_adv_v || do_adv_h) && _advection!(Xp, wzp, fields, ws, timestate)
        do_conv && _convergence!(Xp, fields, timestate, ws)

        @inbounds for j in 1:ydim
            @turbo for i in 1:xdim
                Xp[i+nghost, j] += dX_diff[i, j] + dX_adv[i, j] + dX_conv[i, j]
            end
        end
        refresh_ghosts!(Xp)
    end

    # Final difference
    @inbounds for j in 1:ydim
        @simd for i in 1:xdim
            dX_out[i, j] = Xp[i+nghost, j] - X_in[i, j]
        end
    end

    return nothing
end
