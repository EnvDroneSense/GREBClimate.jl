# ── Grid dimensions ───────────────────────────────────────────────
"Number of longitude grid points."
const xdim = 96
"Number of latitude grid points."
const ydim = 48
const dlon = 360.0f0 / xdim                     # longitude spacing [degrees]
const dlat = 180.0f0 / ydim                     # latitude spacing  [degrees]

# ── Time stepping (the calendar is in calendar.jl) ────────────────
const Δt = 12.0f0 * 3600.0f0                    # main time step [s] (12 hours)
const Δt_crcl = 1800.0f0                        # circulation sub-time step [s] (30 min)
const ntime = max(1, Int(round(Δt / Δt_crcl)))  # Number of sub-steps within one main time step

# ── Numerical limits ──────────────────────────────────────────────
const min_T_K = 40.0f0              # numerical-stability floor [K]
const max_humidity_change = 0.020f0 # Maximum humidity increment [kg/kg]
const min_humidity_change = 0.9f0   # Fraction of humidity that can be removed

# ── Natural constants ────────────────────────────────────────────
const pi_f32 = Float32(pi)   # π (model precision)
const σ = 5.6704f-8            # Stefan-Boltzmann constant [W/m²/K⁴]
const ρ_ocean = 999.1f0        # density of water at T=15°C [kg/m³]
const ρ_land = 2600.0f0        # density of solid rock [kg/m³]
const ρ_air = 1.2f0            # density of air at 20°C at sea
const grav = 9.81f0            # gravitational acceleration [m/s²]
const cp_ocean = 4186.0f0      # specific heat of water at T=15°C [J/kg/K]
const cp_land = cp_ocean / 4.5f0 # specific heat of dry land [J/kg/K]
const cp_air = 1005.0f0        # specific heat of air [J/kg/K]

# ── Column depths [m] ───────────────────────────────────────────
const d_ocean = 50.0f0                      # ocean column
const d_land = 2.0f0                        # land column
const d_air = 5000.0f0                      # air column

# ── Heat capacities [J/K/m²] ────────────────────────────────────
const cap_ocean = cp_ocean * ρ_ocean        # 1m ocean
const cap_land = cp_land * ρ_land * d_land  # land column
const cap_air = cp_air * ρ_air * d_air      # air column
const ΔT_air_factor = Δt / cap_air           # [K m²/W]: flux to air temperature change per step

# ── Sensible heat [W/K/m²] ──────────────────────────────────────
const ct_sens = 22.5f0                       # sensible heat coupling

# ── Albedo parameters ───────────────────────────────────────────
const da_ice = 0.25f0                        # albedo increase for ice-cover
const a_no_ice = 0.1f0                       # albedo ice-free
const a_cloud = 0.35f0                       # cloud albedo

# ── Ice/snow temperature thresholds [K] ─────────────────────────
const Tl_ice1 = 273.15f0 - 10.0f0            # land: full ice albedo
const Tl_ice2 = 273.15f0                     # land: no ice albedo
const To_ice1 = 273.15f0 - 7.0f0             # ocean: full ice
const To_ice2 = 273.15f0 - 1.7f0             # ocean: no ice albedo

# Precomputed inverse ranges (avoids division in hot loops)
const inv_To_ice_range = 1.0f0 / (To_ice2 - To_ice1)
const inv_Tl_ice_range = 1.0f0 / (Tl_ice2 - Tl_ice1)

# Ice fraction at temperature `T`: 1 at or below `T1`, 0 at or above `T2`,
# linear in between. A macro so that the `@turbo` loops see the expression.
macro ice_ramp(T, T1, T2, inv_range)
    return esc(:(ifelse($T <= $T1, 1.0f0, ifelse($T >= $T2, 0.0f0, 1.0f0 - ($T - $T1) * $inv_range))))
end

# ── Deep ocean ──────────────────────────────────────────────────
const co_turb = 5.0f0                        # turbulent mixing coefficient [W/K/m²]
const c_effmix = 0.5f0                       # mixing efficiency
const turb_coeff = Δt * co_turb / cap_ocean  # precomputed mixing coefficient

# ── Atmospheric transport ───────────────────────────────────────
const κ = 8f5                          # diffusion coefficient [m²/s]

# ── Latent heat / hydrology ─────────────────────────────────────
const ce = 2f-3                        # latent heat transfer coefficient
const cq_latent = 2.257f6              # latent heat of evaporation [J/kg]
const cq_rain = -0.1f0 / 24.0f0 / 3600.0f0   # rain-related vapor decrease [1/s]

# ── Scaling heights [m] ─────────────────────────────────────────
const z_air = 8400.0f0                 # heat & CO2 scaling height
const z_vapor = 5000.0f0               # water vapor scaling height
const convergence_factor = Δt_crcl / z_vapor * 2.5f0 / (ρ_air * grav)

# ── Regression factor [kg/m³] ───────────────────────────────────
const r_qviwv = 2.6736f3               # VIWV ↔ q_air regression factor
const q_to_mm_per_day = r_qviwv * 86400.0f0  # kg/kg → mm/day conversion

# ── solar factor [%] ────────────────────────────────────────────
const solar_percent = 100.0f0        # default 100%

# ── Emissivity parameters (LWradiation!) ──────────────────────────
const emissivity_fit = Float32[9.0721, 106.7252, 61.5562, 0.0179, 0.0028,
                      0.0570, 0.3462, 2.3406, 0.7032, 1.0662]

# ── Transport geometry [m, degrees] ───────────────────────────────
const deg_grid = 2.0f0 * pi_f32 * 6.371f6 / 360.0f0
const dyy_grid = dlat * deg_grid
const lat_grid = Float32[dlat * k - dlat / 2.0f0 - 90.0f0 for k in 1:ydim]
const dxlat_grid = Float32[dlon * deg_grid * cos(2.0f0 * pi_f32 / 360.0f0 * lat_grid[k]) for k in 1:ydim]

# ── Diffusion coefficients ──────────────────────────────────────
const ccy_diff = κ * Δt_crcl / dyy_grid^2
const ccx_diff = Float32[κ * Δt_crcl / dxlat_grid[k]^2 for k in 1:ydim]

# ── Advection coefficients ──────────────────────────────────────
const ccy_adv = Δt_crcl / dyy_grid / 2.0f0
const ccx_adv = Float32[Δt_crcl / dxlat_grid[k] / 2.0f0 for k in 1:ydim]

# ── Periodic ghost cells (longitude) ──────────────────────────────
"""
    nghost

*Periodic ghost cells* held at each end of the longitude axis by the
circulation buffers, so wrap-around is a unit-stride read instead of a gather.

A ghosted column stores `A[i, k]` at row `i + nghost`, with each end mirroring
the other:

    row    1 … 3  │  4 … 99         │  100 … 102
    holds  94…96  │  1 … 96 (real)  │  1 … 3

So `A[j-1, k]` is just `P[j+2, k]`.

`nghost` is the zonal stencils' reach, so it is 3 and they spell their offsets
out literally (`j` … `j+6`); widening the stencil means widening both.
"""
const nghost = 3
"Row count of a ghosted circulation buffer (`xdim + 2nghost`)."
const xghost = xdim + 2 * nghost

# More ghost cells than grid points would wrap past the opposite edge.
@assert nghost <= xdim

# ── Polar rows: sub-stepping of diffusion! and advection! ─────────
const polar_threshold = 2.5f5  # 250 km in meters
const is_polar = [dxlat_grid[k] <= polar_threshold for k in 1:ydim]

function _polar_diff_step(k)
    dd = max(1, round(Int, Δt_crcl / (dxlat_grid[k]^2 / κ)))
    dtdff2 = Δt_crcl / dd
    time2 = max(1, round(Int, Δt_crcl / dtdff2))
    return (time2=time2, ccx2=κ * dtdff2 / dxlat_grid[k]^2)
end
function _polar_adv_step(k)
    dd = max(1, round(Int, Δt_crcl / (dxlat_grid[k] / 10.0f0)))
    dtdff2 = Δt_crcl / dd
    time2 = max(1, round(Int, Δt_crcl / dtdff2))
    return (time2=time2, ccx2=dtdff2 / dxlat_grid[k] / 2.0f0)
end

const polar_diff_time2 = [_polar_diff_step(k).time2 for k in 1:ydim]
const polar_diff_ccx2 = Float32[_polar_diff_step(k).ccx2 for k in 1:ydim]
const polar_adv_time2 = [_polar_adv_step(k).time2 for k in 1:ydim]
const polar_adv_ccx2 = Float32[_polar_adv_step(k).ccx2 for k in 1:ydim]
