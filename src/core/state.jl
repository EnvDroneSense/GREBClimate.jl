"""
    CirculationWorkspace

Pre-allocated buffers for diffusion, advection, and circulation calculations.
Reused across all time steps to eliminate allocations.
"""
Base.@kwdef mutable struct CirculationWorkspace
    # Polar sub-stepping buffers. `T1h` carries ghost cells (see `nghost`).
    T1h::Vector{Float32} = zeros(Float32, xghost)  # polar sub-stepping (ghosted)
    dTxh::Vector{Float32} = zeros(Float32, xdim)  # polar increment (Jacobi scratch)

    # Circulation work arrays. `X_work` and `wz_ghost` carry ghost cells
    # (`xghost, ydim`) so the zonal stencils load neighbours contiguously.
    X_work::Matrix{Float32} = zeros(Float32, xghost, ydim)  # circulation work array (ghosted)
    wz_ghost::Matrix{Float32} = zeros(Float32, xghost, ydim)  # ghosted copy of wz_air/wz_vapor
    dX_diff::Matrix{Float32} = zeros(Float32, xdim, ydim)  # diffusion output
    dX_adv::Matrix{Float32} = zeros(Float32, xdim, ydim)  # advection output
    dX_conv::Matrix{Float32} = zeros(Float32, xdim, ydim)  # convection output

    # Tendency buffers
    Q_sens_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Sensible heat flux buffer

    # State buffers
    Ts0_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Surface temperature output
    Ta0_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Air temperature output
    To0_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Ocean temperature output
    q0_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Humidity output

    # LW radiation buffers
    e_co2_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # spatial CO2 buffer
    e_vapor_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # spatial water vapor buffer
    em_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # spatial emissivity buffer
    LW_surf_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Surface longwave
    LW_down_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Downwelling longwave
    LW_up_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Upwelling longwave

    # Hydrology
    Q_lat_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)
    Q_lat_air_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)
    dq_eva_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)
    dq_rain_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)

    # Deep_ocean
    dT_ocean_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)
    dTo_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)

    # Dedicated circulation output
    dTa_crcl::Matrix{Float32} = zeros(Float32, xdim, ydim)  # temperature tendency
    dq_crcl::Matrix{Float32} = zeros(Float32, xdim, ydim)  # humidity tendency

    # SWradiation
    ice_cover_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # ice fraction
    a_surf_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # surface albedo
    albedo_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # combined albedo (surface + atmosphere)
    sw_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # net shortwave flux

    # time_loop
    precip_out::Matrix{Float32} = zeros(Float32, xdim, ydim)  # precipitation output
    evap_out::Matrix{Float32} = zeros(Float32, xdim, ydim)  # evaporation output
    qcrcl_out::Matrix{Float32} = zeros(Float32, xdim, ydim)  # circulation moisture output
    term_north::Vector{Float32} = zeros(Float32, xdim)  # diffusion term of the northernmost row
    term_south::Vector{Float32} = zeros(Float32, xdim)  # diffusion term of the southernmost row
end

"""
    SurfaceState

A run's current surface/atmosphere state - `Ts`, `Ta`, `To`, `q` - passed as
one argument to [`diagnostics!`](@ref), [`output!`](@ref), [`time_loop!`](@ref),
and [`qflux_correction!`](@ref). A thin reference wrapper around
already-allocated arrays; construct once per run/call (like `ws`/`acc`),
never inside the per-timestep loop.
"""
struct SurfaceState
    Ts::Matrix{Float32}
    Ta::Matrix{Float32}
    To::Matrix{Float32}
    q::Matrix{Float32}
end

"""
    MonthlyAccumulator

Accumulates fields over a month for monthly-mean output. Reset after each
month via `reset!`. `olrmm` and `lwdownmm` are summed with the signs of
[`MonthlyRecord`](@ref) (`olr` positive upward, `lwdown` positive into the
surface).
"""
Base.@kwdef mutable struct MonthlyAccumulator
    Tmm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Surface temperature accumulator
    Tamm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Air temperature accumulator
    Tomm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Ocean temperature accumulator
    qmm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Humidity accumulator
    apmm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Albedo accumulator
    icemm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Ice fraction accumulator
    precipmm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Precipitation accumulator
    evapmm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Evaporation accumulator
    qcrclmm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Circulation moisture accumulator
    swmm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Shortwave radiation accumulator
    lwmm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Longwave radiation accumulator
    qlatmm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Latent heat accumulator
    qsensmm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Sensible heat accumulator
    olrmm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Outgoing longwave accumulator
    lwdownmm::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Downward longwave accumulator
end

function reset!(acc::MonthlyAccumulator)
    fill!(acc.Tmm, 0.0f0)
    fill!(acc.Tamm, 0.0f0)
    fill!(acc.Tomm, 0.0f0)
    fill!(acc.qmm, 0.0f0)
    fill!(acc.apmm, 0.0f0)
    fill!(acc.icemm, 0.0f0)
    fill!(acc.precipmm, 0.0f0)
    fill!(acc.evapmm, 0.0f0)
    fill!(acc.qcrclmm, 0.0f0)
    fill!(acc.swmm, 0.0f0)
    fill!(acc.lwmm, 0.0f0)
    fill!(acc.qlatmm, 0.0f0)
    fill!(acc.qsensmm, 0.0f0)
    fill!(acc.olrmm, 0.0f0)
    fill!(acc.lwdownmm, 0.0f0)
end

function accumulate!(acc::MonthlyAccumulator, Ts, Ta, To, q, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens,
    lw_up, lw_down, em)
    Tmm = acc.Tmm; Tamm = acc.Tamm; Tomm = acc.Tomm; qmm = acc.qmm
    apmm = acc.apmm; icemm = acc.icemm
    precipmm = acc.precipmm; evapmm = acc.evapmm; qcrclmm = acc.qcrclmm
    swmm = acc.swmm; lwmm = acc.lwmm; qlatmm = acc.qlatmm; qsensmm = acc.qsensmm
    olrmm = acc.olrmm; lwdownmm = acc.lwdownmm

    @turbo for j in 1:ydim
        for i in 1:xdim
            Tmm[i, j] += Ts[i, j]
            Tamm[i, j] += Ta[i, j]
            Tomm[i, j] += To[i, j]
            qmm[i, j] += q[i, j]
            apmm[i, j] += albedo[i, j]
            icemm[i, j] += ice[i, j]
            precipmm[i, j] += precip[i, j]
            evapmm[i, j] += evap[i, j]
            qcrclmm[i, j] += qcrcl[i, j]
            swmm[i, j] += sw[i, j]
            lwmm[i, j] += lw[i, j]
            qlatmm[i, j] += qlat[i, j]
            qsensmm[i, j] += qsens[i, j]
            olrmm[i, j] -= lw_up[i, j] + (1.0f0 - em[i, j]) * lw[i, j]
            lwdownmm[i, j] -= lw_down[i, j]
        end
    end
end

"""
    ClimateFields

Loaded climatology, derived grid fields, flux corrections, the CO2 mask and
the insolation table: what `load_greb_jld2!` fills in and every physics
function reads. It is passed as an argument, never held as global state. One
loaded instance can serve several runs one after another: [`greb_model!`](@ref)
restores the input fields it changes and derives the others again at the
start of each run.

`ClimateFields()` builds an all-zero instance with `loaded = false`;
`load_greb_jld2!` sets `loaded = true`. [`greb_model!`](@ref) refuses unloaded
fields unless `allow_uninitialized=true`: an all-zero climatology raises no
error, it runs and returns NaN in every output field.
"""
Base.@kwdef mutable struct ClimateFields
    # 2D fields (xdim, ydim)
    z_topo::Matrix{Float32} = zeros(Float32, xdim, ydim)  # topography [m] (<0: ocean)
    glacier::Matrix{Float32} = zeros(Float32, xdim, ydim)  # glacier mask (>0.5: glacier)
    z_ocean::Matrix{Float32} = zeros(Float32, xdim, ydim)  # derived ocean depth [m]
    cap_surf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # surface heat capacity [J/K/m²]
    wz_air::Matrix{Float32} = zeros(Float32, xdim, ydim)  # exp(-z_topo / z_air)
    wz_vapor::Matrix{Float32} = zeros(Float32, xdim, ydim)  # exp(-z_topo / z_vapor)
    rain_limit::Matrix{Float32} = zeros(Float32, xdim, ydim)  # -0.0015/(wz_vapor*r_qviwv*86400)
    # 3D climate fields (xdim, ydim, nstep_yr)
    Tclim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # surface temperature [K]
    uclim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # zonal wind [m/s]
    vclim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # meridional wind [m/s]
    qclim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # atmospheric humidity [kg/kg]
    mldclim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # mixed-layer depth [m]
    omegaclim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # vertical velocity [Pa/s]
    omegastdclim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # omega std deviation [Pa/s]
    wsclim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # wind speed [m/s]

    # Anomaly fields for ENSO/climate-change experiments
    Tclim_anom_enso::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    uclim_anom_enso::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    vclim_anom_enso::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    omegaclim_anom_enso::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    wsclim_anom_enso::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    Tclim_anom_cc::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    uclim_anom_cc::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    vclim_anom_cc::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    omegaclim_anom_cc::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    wsclim_anom_cc::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)

    # The winds split by sign (derive_fields!): pos + neg = the wind
    uclim_pos::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # eastward part, 0 elsewhere
    uclim_neg::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # westward part, 0 elsewhere
    vclim_pos::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # northward part, 0 elsewhere
    vclim_neg::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # southward part, 0 elsewhere

    Toclim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # deep ocean temperature [K]
    cldclim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # cloud cover fraction
    swetclim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # soil wetness [0-1]

    # Solar / radiation
    sw_solar::Matrix{Float32} = zeros(Float32, ydim, nstep_yr)  # 24hr mean solar radiation [W/m²] (ydim, nstep_yr)
    dTrad::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # Tatmos-radiation offset

    # Flux correction arrays (zeros unless loaded from file)
    TF_correct::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    qF_correct::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    ToF_correct::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)

    # Regional CO2 mask (1.0 = full CO2, 0.5 = half CO2)
    co2_part::Matrix{Float32} = ones(Float32, xdim, ydim)

    # false for a bare `ClimateFields()`; set by `load_greb_jld2!`. See the
    # docstring above and `greb_model!`'s `allow_uninitialized` keyword.
    loaded::Bool = false
end

"""
    TimeState

Tracks the model's position within the current year: `jday` (calendar
day, 1..365) and `ityr` (timestep-of-year, 1..`nstep_yr`). Mutated in
place each timestep by [`time_loop!`](@ref)/[`qflux_correction!`](@ref).
"""
mutable struct TimeState
    jday::Int  # Current calendar day in year [1..365]
    ityr::Int  # Current timestep in year [1..730]
end

# The land test as a macro for the `@turbo` loops, which compile a comparison
# they can see better than a function call
macro is_land(z)
    return esc(:($z > 0.0f0))
end

"""
    is_land(z_topo)

Whether a cell of topographic height `z_topo` [m] is land: above 0 m.
Everything else is ocean. The one land/ocean test of the model.
"""
is_land(z_topo) = @is_land(z_topo)

"""
    derive_fields!(fields::ClimateFields, p::Processes)

Compute every field of `fields` that follows from its input maps: the
radiation-temperature offset `dTrad` (from `Tclim`), the deep-ocean depth
`z_ocean` (from `mldclim`), the pressure weights `wz_air`/`wz_vapor` and the
rain limit (from `z_topo`), the surface heat capacity `cap_surf` (from
`z_topo` and `mldclim`) and the winds split by sign. Call it again after
changing an input map; [`init_model!`](@ref) calls it once per run.
"""
function derive_fields!(fields::ClimateFields, p::Processes)
    z_topo = fields.z_topo
    mldclim = fields.mldclim

    @. fields.dTrad = -0.16f0 * fields.Tclim - 5.0f0
    # Three times the deepest mixed layer of the year
    fields.z_ocean .= 3.0f0 .* dropdims(maximum(mldclim; dims=3); dims=3)

    @. fields.wz_air = exp(-z_topo / z_air)
    @. fields.wz_vapor = exp(-z_topo / z_vapor)
    @. fields.rain_limit = -0.0015f0 / (fields.wz_vapor * r_qviwv * 86400.0f0)

    cap_surf = fields.cap_surf
    for j in 1:ydim, i in 1:xdim
        cap_surf[i, j] = is_land(z_topo[i, j]) || p.ocean === :none ? cap_land : cap_ocean * mldclim[i, j, 1]
    end

    split_winds!(fields)
    return fields
end

# Upwind advection reads the eastward/westward and northward/southward parts
function split_winds!(fields::ClimateFields)
    _split_sign!(fields.uclim_pos, fields.uclim_neg, fields.uclim)
    _split_sign!(fields.vclim_pos, fields.vclim_neg, fields.vclim)
    return fields
end

function _split_sign!(pos, neg, x)
    @turbo for i in eachindex(x)
        pos[i] = ifelse(x[i] >= 0.0f0, x[i], 0.0f0)
        neg[i] = ifelse(x[i] < 0.0f0, x[i], 0.0f0)
    end
    return nothing
end

"""
    ModelState

Per-run mutable state that isn't climatology: the runtime solar-forcing
multiplier (`SWradiation!` reads it) and the surface-temperature accumulator
behind the annual progress line (`diagnostics!` reads/writes it). One instance
per `greb_model!` run.

This is scratch space for the printed summary, not an output path - `Tsmn` is
averaged, printed and zeroed within a single `diagnostics!` call, so it never
holds a readable annual mean once the call returns. Model output is the
`Vector{MonthlyRecord}` that [`greb_model!`](@ref) returns.
"""
mutable struct ModelState
    sw_solar_forcing::Float32   # runtime solar multiplier used by SWradiation!
    Tsmn::Matrix{Float32}       # surface-temperature accumulator for the progress line
end

function ModelState()
    ModelState(1.0f0, zeros(Float32, xdim, ydim))
end

"""
    MonthlyRecord

One monthly-mean output record: a `NamedTuple` with fields `Ts`, `Ta`, `To`,
`q`, `albedo`, `ice`, `precip`, `evap`, `qcrcl`, `sw`, `lw`, `qlat`, `qsens`,
`olr`, `lwdown`, each an `(xdim, ydim)` `Matrix{Float32}`. Produced by
[`output!`](@ref); `greb_model!`'s `ctrl`/`scnr` results are
`Vector{MonthlyRecord}`.

The fluxes are in W/m2. `sw`, `lw`, `qlat` and `qsens` are positive into the
surface; `lw` is the surface's own emission alone, so it is negative. `olr` is
the longwave leaving to space, positive upward; `lwdown` is the longwave the
air sends down, positive into the surface.
"""
const MonthlyRecord = NamedTuple{(:Ts, :Ta, :To, :q, :albedo, :ice, :precip, :evap, :qcrcl, :sw, :lw, :qlat, :qsens, :olr, :lwdown),NTuple{15,Matrix{Float32}}};
