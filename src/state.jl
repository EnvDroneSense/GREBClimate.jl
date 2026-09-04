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
    temp_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # general workspace (humidity-update scratch)
    Q_sens_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Sensible heat flux buffer
    crcl::Matrix{Float32} = zeros(Float32, xdim, ydim)  # dq_crcl (zero stand-in when log_crcl_dmc is off)

    # State buffers
    Ts0_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Surface temperature output
    Ta0_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Air temperature output
    To0_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Ocean temperature output
    q0_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Humidity output

    # LW radiation buffers
    e_co2_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # spatial CO₂ buffer
    e_vapor_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # spatial water vapor buffer
    em_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # spatial emissivity buffer
    LW_surf_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Surface longwave
    LW_down_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Downwelling longwave
    LW_up_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Upwelling longwave

    # Hydrology buffers
    qs::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Saturation humidity buffer
    Tskin::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Skin temperature buffer
    rq::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Relative humidity buffer
    ws_base::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Base wind speed buffer
    # Hydrology
    Q_lat_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)
    Q_lat_air_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)
    dq_eva_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)
    dq_rain_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)
    cE_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Surface exchange coefficient buffer

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
    a_atmos_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # atmospheric albedo
    sw_buf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # net shortwave flux

    # time_loop
    precip_out::Matrix{Float32} = zeros(Float32, xdim, ydim)  # precipitation output
    evap_out::Matrix{Float32} = zeros(Float32, xdim, ydim)  # evaporation output
    qcrcl_out::Matrix{Float32} = zeros(Float32, xdim, ydim)  # circulation moisture output
    term_north::Vector{Float32} = zeros(Float32, xdim)  # northern boundary term
    term_south::Vector{Float32} = zeros(Float32, xdim)  # southern boundary term
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

Accumulates fields over a month for monthly-mean output.
Reset after each month via `reset!`.
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
end

function accumulate!(acc::MonthlyAccumulator, Ts, Ta, To, q, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens)
    Tmm = acc.Tmm; Tamm = acc.Tamm; Tomm = acc.Tomm; qmm = acc.qmm
    apmm = acc.apmm; icemm = acc.icemm
    precipmm = acc.precipmm; evapmm = acc.evapmm; qcrclmm = acc.qcrclmm
    swmm = acc.swmm; lwmm = acc.lwmm; qlatmm = acc.qlatmm; qsensmm = acc.qsensmm

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
        end
    end
end

"""
    ClimateFields

Loaded climatology, derived grid fields, flux corrections, and the
regional-CO2 mask/solar table - everything `load_greb_jld2!` fills in and
every physics function reads. One instance per `greb_model!` run; never
shared as global state.

`ClimateFields()` builds an all-zero instance and leaves `loaded = false`;
`load_greb_jld2!` sets `loaded = true` once real climatology is in place.
[`greb_model!`](@ref) refuses to run unloaded fields unless explicitly
told to via `allow_uninitialized=true` - an all-zero climatology produces
a physically meaningless world pinned at the 40 K floor (~-233 °C) rather
than an error, so the
flag exists to keep that path opt-in.
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

    # Precomputed wind sign splits
    uclim_m::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # negative u components
    uclim_p::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # positive u components
    vclim_m::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # negative v components
    vclim_p::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # positive v components

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

    # Regional CO₂ mask (1.0 = full CO₂, 0.5 = half CO₂)
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
each an `(xdim, ydim)` `Matrix{Float32}`. Produced by [`output!`](@ref);
`greb_model!`'s `ctrl`/`scnr` results are `Vector{MonthlyRecord}`.
"""
const MonthlyRecord = NamedTuple{(:Ts, :Ta, :To, :q, :albedo, :ice, :precip, :evap, :qcrcl, :sw, :lw, :qlat, :qsens),NTuple{13,Matrix{Float32}}};
