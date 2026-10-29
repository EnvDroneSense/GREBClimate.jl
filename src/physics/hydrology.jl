const _qsat_scale = 3.75f-3
const _qsat_a = 17.08085f0
const _qsat_b = 234.175f0
const _hydro_gust_land = 4.0f0
const _hydro_gust_ocean = 9.0f0
const _hydro_ce_land = 0.25f0 * ce
const _hydro_ce_ocean = 0.58f0 * ce
const _hydro_latent_factor = cq_latent * ρ_air * ce
# Squared gust speeds [m²/s²] and drag coefficients of the other three schemes
const _skin_gust_land = 132.25f0
const _skin_gust_ocean = 29.16f0
const _orig_gust_land = _hydro_gust_land + 144.0f0
const _orig_gust_ocean = _hydro_gust_ocean + 50.41f0  # 7.1^2
const _orig_coeff_land = 0.04f0
const _orig_coeff_ocean = 0.73f0
const _skin_gust2_land = 81.0f0  # 9.0^2
const _skin_gust2_ocean = 16.0f0  # 4.0^2
const _skin_coeff_land = 0.56f0
const _skin_coeff_ocean = 0.79f0

"""
    hydro!(ws::ModelWorkspace, Ts, q, fields::ClimateFields, slice, p::Processes, h::ResolvedHydrology)

Compute latent heat flux and evaporation/rain tendencies. `h.evaporation`
selects the evaporation scheme; `h.rain` and its coefficients
(`c_q`/`c_rq`/`c_omega`/`c_omega_std`) the rain regression. Returns zeros
without an atmosphere or with `p.hydrology` other than `:full`.
Returns `(Q_lat, Q_lat_air, dq_eva, dq_rain)`. The saturation formula is not
finite for `Ts` at or below 38.975 K; a run stays above it through the 40 K
floor.
"""
function hydro!(ws::ModelWorkspace, Ts, q, fields::ClimateFields, slice::Int, p::Processes, h::ResolvedHydrology)
    c_q = h.c_q
    c_rq = h.c_rq
    c_omega = h.c_omega
    c_omega_std = h.c_omega_std

    if !p.atmosphere || p.hydrology !== :full
        fill!(ws.Q_lat, 0.0f0)
        fill!(ws.Q_lat_air, 0.0f0)
        fill!(ws.dq_eva, 0.0f0)
        fill!(ws.dq_rain, 0.0f0)
        return (Q_lat=ws.Q_lat, Q_lat_air=ws.Q_lat_air,
            dq_eva=ws.dq_eva, dq_rain=ws.dq_rain)
    end

    z_topo = fields.z_topo
    wz_air = fields.wz_air
    u = clim_slice(fields.u_clim, slice)
    v = clim_slice(fields.v_clim, slice)
    swet = clim_slice(fields.soil_wetness_clim, slice)
    omega = clim_slice(fields.omega_clim, slice)
    omega_std = clim_slice(fields.omega_std_clim, slice)
    rain_limit = fields.rain_limit
    apply_rain_limit = h.rain === :rh

    Q_lat = ws.Q_lat
    Q_lat_air = ws.Q_lat_air
    dq_eva = ws.dq_eva
    dq_rain = ws.dq_rain

    # Saturation humidity, relative humidity, evaporation, precipitation, the
    # optional rain-limit clamp, and water-vapor tendencies
    if h.evaporation === :original
        for j in 1:ydim
            @inbounds @fastmath @simd ivdep for i in 1:xdim
                T = Ts[i, j] - 273.15f0
                qs = max(_qsat_scale * fastexp(_qsat_a * T / (T + _qsat_b)) * wz_air[i, j], 1.0f-8)
                rq = q[i, j] / qs

                u_val = u[i, j]
                v_val = v[i, j]
                wind = sqrt(u_val*u_val + v_val*v_val)
                wind = sqrt(wind*wind + ifelse(@is_land(z_topo[i, j]), _hydro_gust_land, _hydro_gust_ocean))
                qlat = (q[i, j] - qs) * wind * _hydro_latent_factor * swet[i, j]
                Q_lat[i, j] = qlat

                drain = (c_q + c_rq * rq + c_omega * omega[i, j] + c_omega_std * omega_std[i, j]) * cq_rain * q[i, j]
                limit_val = rain_limit[i, j]
                drain = ifelse(apply_rain_limit & (drain >= limit_val), limit_val, drain)
                dq_rain[i, j] = drain

                dq_eva[i, j] = -qlat / cq_latent / r_qviwv
                Q_lat_air[i, j] = -drain * cq_latent * r_qviwv
            end
        end
    elseif h.evaporation === :skin
        wind_speed = clim_slice(fields.wind_speed_clim, slice)
        for j in 1:ydim
            @inbounds @fastmath @simd ivdep for i in 1:xdim
                T0 = Ts[i, j] - 273.15f0
                qs0 = max(_qsat_scale * fastexp(_qsat_a * T0 / (T0 + _qsat_b)) * wz_air[i, j], 1.0f-8)
                rq = q[i, j] / qs0

                Tskin = ifelse(@is_land(z_topo[i, j]), Ts[i, j] + 5.0f0, Ts[i, j] + 1.0f0)
                Tskin = ifelse(Tskin < 200.0f0, 200.0f0, Tskin)
                T = Tskin - 273.15f0
                qs_val = _qsat_scale * fastexp(_qsat_a * T / (T + _qsat_b)) * wz_air[i, j]

                mean_wind = wind_speed[i, j]
                gust = ifelse(@is_land(z_topo[i, j]), _skin_gust_land, _skin_gust_ocean)
                wind = sqrt(mean_wind*mean_wind + gust)

                cE = ifelse(@is_land(z_topo[i, j]), _hydro_ce_land, _hydro_ce_ocean)
                qlat = cE * wind * ρ_air * cq_latent * (q[i, j] - qs_val) * swet[i, j]
                Q_lat[i, j] = qlat

                drain = (c_q + c_rq * rq + c_omega * omega[i, j] + c_omega_std * omega_std[i, j]) * cq_rain * q[i, j]
                limit_val = rain_limit[i, j]
                drain = ifelse(apply_rain_limit & (drain >= limit_val), limit_val, drain)
                dq_rain[i, j] = drain

                dq_eva[i, j] = -qlat / cq_latent / r_qviwv
                Q_lat_air[i, j] = -drain * cq_latent * r_qviwv
            end
        end
    elseif h.evaporation === :original_gust
        for j in 1:ydim
            @inbounds @fastmath @simd ivdep for i in 1:xdim
                T = Ts[i, j] - 273.15f0
                qs = max(_qsat_scale * fastexp(_qsat_a * T / (T + _qsat_b)) * wz_air[i, j], 1.0f-8)
                rq = q[i, j] / qs

                u_val = u[i, j]
                v_val = v[i, j]
                wind = sqrt(u_val*u_val + v_val*v_val)
                wind = sqrt(wind*wind + ifelse(@is_land(z_topo[i, j]), _orig_gust_land, _orig_gust_ocean))
                coeff = ifelse(@is_land(z_topo[i, j]), _orig_coeff_land, _orig_coeff_ocean)
                qlat = (q[i, j] - qs) * wind * cq_latent * ρ_air * coeff * ce * swet[i, j]
                Q_lat[i, j] = qlat

                drain = (c_q + c_rq * rq + c_omega * omega[i, j] + c_omega_std * omega_std[i, j]) * cq_rain * q[i, j]
                limit_val = rain_limit[i, j]
                drain = ifelse(apply_rain_limit & (drain >= limit_val), limit_val, drain)
                dq_rain[i, j] = drain

                dq_eva[i, j] = -qlat / cq_latent / r_qviwv
                Q_lat_air[i, j] = -drain * cq_latent * r_qviwv
            end
        end
    elseif h.evaporation === :skin_gust
        wind_speed = clim_slice(fields.wind_speed_clim, slice)
        for j in 1:ydim
            @inbounds @fastmath @simd ivdep for i in 1:xdim
                T = Ts[i, j] - 273.15f0
                qs = max(_qsat_scale * fastexp(_qsat_a * T / (T + _qsat_b)) * wz_air[i, j], 1.0f-8)
                rq = q[i, j] / qs

                wind = wind_speed[i, j]
                wind = sqrt(wind*wind + ifelse(@is_land(z_topo[i, j]), _skin_gust2_land, _skin_gust2_ocean))
                coeff = ifelse(@is_land(z_topo[i, j]), _skin_coeff_land, _skin_coeff_ocean)
                qlat = (q[i, j] - qs) * wind * cq_latent * ρ_air * coeff * ce * swet[i, j]
                Q_lat[i, j] = qlat

                drain = (c_q + c_rq * rq + c_omega * omega[i, j] + c_omega_std * omega_std[i, j]) * cq_rain * q[i, j]
                limit_val = rain_limit[i, j]
                drain = ifelse(apply_rain_limit & (drain >= limit_val), limit_val, drain)
                dq_rain[i, j] = drain

                dq_eva[i, j] = -qlat / cq_latent / r_qviwv
                Q_lat_air[i, j] = -drain * cq_latent * r_qviwv
            end
        end
    else
        error("Unknown evaporation scheme :$(h.evaporation)")
    end

    return (Q_lat=Q_lat,
        Q_lat_air=Q_lat_air,
        dq_eva=dq_eva,
        dq_rain=dq_rain)
end
