"""
    seaice!(Ts0, fields::ClimateFields, slice, p::Processes)

Updates `fields.cap_surf` (surface heat capacity) for ocean points based on
`Ts0`-derived ice fraction, blending land/open-ocean/ice capacities. No-op
with `p.ocean = :none`; skips the ice blend with `p.ice_albedo = false`.
"""
function seaice!(Ts0, fields::ClimateFields, slice::Int, p::Processes)
    mld = clim_slice(fields.mld_clim, slice)
    z_topo = fields.z_topo
    glacier = fields.glacier
    cap_surf = fields.cap_surf

    if p.ocean === :none
        return nothing   # No ice feedback: skip sea ice calculation
    end

    if p.ice_albedo
        # Compute ice‑dependent heat capacity for ocean points
        @turbo for i in 1:xdim, j in 1:ydim
            is_ocean = !@is_land(z_topo[i, j])
            T = Ts0[i, j]
            mld_val = mld[i, j]
            cap_open = cap_ocean * mld_val

            # Ice fraction (0 = no ice, 1 = full ice)
            ice_frac = @ice_ramp(T, To_ice1, To_ice2, inv_To_ice_range)

            # Blend between land (ice) and open ocean capacities
            cap_with_ice = cap_land * ice_frac + cap_open * (1.0f0 - ice_frac)

            # Apply only to ocean points; keep land points unchanged
            cap_surf[i, j] = ifelse(is_ocean, cap_with_ice, cap_surf[i, j])
        end
    else
        # Without the ice-albedo feedback the ocean keeps its open-water capacity
        @turbo for i in 1:xdim, j in 1:ydim
            cap_surf[i, j] = ifelse(@is_land(z_topo[i, j]), cap_land, cap_ocean * mld[i, j])
        end
    end

    # Glacier override: ice sheets have land heat capacity.
    @. cap_surf = ifelse(glacier > 0.5f0, cap_land, cap_surf)
    return nothing
end

"""
    deep_ocean!(Ts, To, fields::ClimateFields, slice, p::Processes, ws::ModelWorkspace)

Computes surface/deep-ocean coupling tendencies (`dT_ocean`, `dTo`) from
mixed-layer-depth entrainment/detrainment and turbulent mixing, active only
where the point is ocean and above the sea-ice threshold. Returns zeros
unless `p.ocean` is `:full`.
"""
function deep_ocean!(Ts, To, fields::ClimateFields, slice::Int, p::Processes, ws::ModelWorkspace)
    # Use pre-allocated zero buffers
    dT_ocean = ws.dT_ocean
    dTo = ws.dTo

    # no deep-ocean coupling
    if p.ocean !== :full
        fill!(dT_ocean, 0.0f0)
        fill!(dTo, 0.0f0)
        return (dT_ocean=dT_ocean, dTo=dTo)
    end

    z_topo = fields.z_topo
    z_ocean = fields.z_ocean

    # ── Change in mixed-layer depth ─────────────────────────
    mld_now = clim_slice(fields.mld_clim, slice)
    mld_prev = clim_slice(fields.mld_clim, previous_slice(slice))

    # Zero buffers first
    fill!(dT_ocean, 0.0f0)
    fill!(dTo, 0.0f0)

    # ── Entrainment & detrainment & turbulent mixing ──────
    @turbo for i in 1:xdim, j in 1:ydim
        is_ocean = !@is_land(z_topo[i, j])
        # Entrainment/detrainment require Ts >= To_ice2
        active = is_ocean & (Ts[i, j] >= To_ice2)
        h_now = mld_now[i, j]
        h_prev = mld_prev[i, j]
        dh = h_now - h_prev
        z_deep = z_ocean[i, j]
        z_rem = z_deep - h_now

        # Entrainment/detrainment contributions (only when active)
        dTo_entr = ifelse(active & (dh < 0.0f0), c_effmix * (-dh / z_rem) *
                                               (Ts[i, j] - To[i, j]), 0.0f0)
        dT_ocean_entr = ifelse(active & (dh > 0.0f0), c_effmix * (dh / h_now) *
                                                    (To[i, j] - Ts[i, j]), 0.0f0)

        # Turbulent mixing (ocean points only, regardless of ice threshold)
        Tx = ifelse(Ts[i, j] > To_ice2, Ts[i, j], To_ice2)
        dTo_turb = ifelse(is_ocean, turb_coeff * (Tx - To[i, j]) / z_rem, 0.0f0)
        dT_ocean_turb = ifelse(is_ocean, turb_coeff * (To[i, j] - Tx) / h_now, 0.0f0)

        # Combine (buffer was zeroed before loop)
        dTo[i, j] = dTo_entr + dTo_turb
        dT_ocean[i, j] = dT_ocean_entr + dT_ocean_turb
    end
    return (dT_ocean=dT_ocean, dTo=dTo)
end
