"""
    SWradiation!(Ts, fields::ClimateFields, state::ModelState, timestate, p::Processes, ws::ModelWorkspace)

Computes ice cover, surface/atmospheric/combined albedo, and net shortwave
flux from `Ts` and the current cloud climatology. Returns
`(SW, albedo, ice_cover)`.
"""
function SWradiation!(Ts, fields::ClimateFields, state::ModelState, timestate, p::Processes, ws::ModelWorkspace)
    # Reuse workspace buffers
    ice_cover = ws.ice_cover # output: ice fraction
    a_surf = ws.a_surf       # surface albedo
    albedo = ws.albedo       # output: combined albedo (surface + atmosphere)
    sw = ws.sw               # output: net shortwave flux

    z_topo = fields.z_topo
    glacier = fields.glacier

    # 1. Ice cover fraction
    @turbo for i in 1:xdim, j in 1:ydim
        T = Ts[i, j]
        land_ice = @ice_ramp(T, Tl_ice1, Tl_ice2, inv_Tl_ice_range)
        ocean_ice = @ice_ramp(T, To_ice1, To_ice2, inv_To_ice_range)
        ice_cover[i, j] = ifelse(@is_land(z_topo[i, j]), land_ice, ocean_ice)
    end

    # 2. Atmospheric albedo
    cld = fields.cloud_clim
    ityr = timestate.ityr

    # 3. Surface albedo
    if p.ice_albedo
        @turbo for i in 1:xdim, j in 1:ydim
            T = Ts[i, j]
            # Ice-free albedo plus the ice increment times the ice fraction
            land_alb = a_no_ice + da_ice * @ice_ramp(T, Tl_ice1, Tl_ice2, inv_Tl_ice_range)
            ocean_alb = a_no_ice + da_ice * @ice_ramp(T, To_ice1, To_ice2, inv_To_ice_range)
            # Choose based on topography
            a_surf[i, j] = ifelse(@is_land(z_topo[i, j]), land_alb, ocean_alb)
            # Glacier override: if glacier mask > 0.5, set to ice albedo
            a_surf[i, j] = ifelse(glacier[i, j] > 0.5f0, a_no_ice + da_ice, a_surf[i, j])
        end
    else
        @. a_surf = a_no_ice
    end

    # 4. albedo + shortwave flux.
    sw_solar = fields.sw_solar
    multiplier = state.sw_solar_forcing * 0.01f0 * solar_percent
    @turbo for j in 1:ydim
        sf = sw_solar[j, ityr] * multiplier
        for i in 1:xdim
            aa = cld[i, j, ityr] * a_cloud
            alb = a_surf[i, j] + aa - a_surf[i, j] * aa
            albedo[i, j] = alb
            sw[i, j] = sf * (1.0f0 - alb)
        end
    end

    return (SW=sw, albedo=albedo, ice_cover=ice_cover)
end

"""
    LWradiation!(Ts, Ta, q, CO2, fields::ClimateFields, timestate, p::Processes, ws::ModelWorkspace)

Computes atmospheric emissivity from CO2/water-vapor/cloud columns, then
surface/upward/downward longwave flux. Without an atmosphere (`p.atmosphere = false`) only
`LW_down` is zeroed - `LW_up` is snapshotted beforehand and keeps its full
value (decouples surface from atmospheric downwelling feedback without
touching the atmosphere's own emission term). Returns
`(LW_surf, LW_up, LW_down, em)`.
"""
function LWradiation!(Ts, Ta, q, CO2, fields::ClimateFields, timestate, p::Processes, ws::ModelWorkspace)
    # Extract workspace buffers
    e_co2 = ws.e_co2      # CO2 [ppm scaled by pressure]
    e_vapor = ws.e_vapor  # water vapor [kg/m²]
    em = ws.em            # emissivity ε_atmos
    LW_surf = ws.LW_surf  # surface long-wave flux [W/m²]
    LW_down = ws.LW_down  # downward long-wave flux [W/m²]
    LW_up = ws.LW_up      # upward long-wave flux [W/m²]

    wz_air = fields.wz_air
    co2_part = fields.co2_part
    ityr = timestate.ityr
    cloud_clim = fields.cloud_clim
    dTrad = fields.dTrad
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = emissivity_fit

    # ── Effective columns, emissivity (log-regression, 10 parameters),
    # cloud adjustment, and surface/downward longwave flux
    @turbo for j in 1:ydim
        for i in 1:xdim
            e_vapor[i, j] = wz_air[i, j] * r_qviwv * q[i, j]
            e_co2[i, j] = wz_air[i, j] * CO2 * co2_part[i, j]
            em_val = p4 * log(p1 * e_co2[i, j] + p2 * e_vapor[i, j] + p3) +
                     p7 +
                     p5 * log(p1 * e_co2[i, j] + p3) +
                     p6 * log(p2 * e_vapor[i, j] + p3)
            em_val = (p8 - cloud_clim[i, j, ityr]) / p9 * (em_val - p10) + p10
            em[i, j] = em_val
            LW_surf[i, j] = -σ * Ts[i, j]^4
            LW_down_val = -em_val * σ * (Ta[i, j] + dTrad[i, j, ityr])^4
            LW_down[i, j] = LW_down_val
            LW_up[i, j] = LW_down_val
        end
    end

    if !p.atmosphere
        LW_down .= 0.0f0
    end

    return (LW_surf=LW_surf, LW_up=LW_up, LW_down=LW_down, em=em)
end
