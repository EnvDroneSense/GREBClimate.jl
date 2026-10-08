"""
    ModelWorkspace

Pre-allocated buffers for diffusion, advection, and circulation calculations.
Reused across all time steps to eliminate allocations.
"""
Base.@kwdef struct ModelWorkspace
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
    Q_sens::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Sensible heat flux buffer

    # State buffers
    Ts0::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Surface temperature output
    Ta0::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Air temperature output
    To0::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Ocean temperature output
    q0::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Humidity output

    # LW radiation buffers
    e_co2::Matrix{Float32} = zeros(Float32, xdim, ydim)  # spatial CO2 buffer
    e_vapor::Matrix{Float32} = zeros(Float32, xdim, ydim)  # spatial water vapor buffer
    em::Matrix{Float32} = zeros(Float32, xdim, ydim)  # spatial emissivity buffer
    LW_surf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Surface longwave
    LW_down::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Downwelling longwave
    LW_up::Matrix{Float32} = zeros(Float32, xdim, ydim)  # Upwelling longwave

    # Hydrology
    Q_lat::Matrix{Float32} = zeros(Float32, xdim, ydim)
    Q_lat_air::Matrix{Float32} = zeros(Float32, xdim, ydim)
    dq_eva::Matrix{Float32} = zeros(Float32, xdim, ydim)
    dq_rain::Matrix{Float32} = zeros(Float32, xdim, ydim)

    # Deep_ocean
    dT_ocean::Matrix{Float32} = zeros(Float32, xdim, ydim)
    dTo::Matrix{Float32} = zeros(Float32, xdim, ydim)

    # Dedicated circulation output
    dTa_crcl::Matrix{Float32} = zeros(Float32, xdim, ydim)  # temperature tendency
    dq_crcl::Matrix{Float32} = zeros(Float32, xdim, ydim)  # humidity tendency

    # SWradiation
    ice_cover::Matrix{Float32} = zeros(Float32, xdim, ydim)  # ice fraction
    a_surf::Matrix{Float32} = zeros(Float32, xdim, ydim)  # surface albedo
    albedo::Matrix{Float32} = zeros(Float32, xdim, ydim)  # combined albedo (surface + atmosphere)
    sw::Matrix{Float32} = zeros(Float32, xdim, ydim)  # net shortwave flux

    # time_loop
    precip::Matrix{Float32} = zeros(Float32, xdim, ydim)  # precipitation output
    evap::Matrix{Float32} = zeros(Float32, xdim, ydim)  # evaporation output
    qcrcl::Matrix{Float32} = zeros(Float32, xdim, ydim)  # circulation moisture output
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
    @output_fields (args...) begin
        name += expr
        name -= expr
    end

Declares the monthly output fields, one row per field. A row says what one
timestep adds to the monthly sum of `name`; each `struct.field` in `expr` is an
`(xdim, ydim)` array, read cell by cell. `args` are the arguments of
`accumulate!` after the accumulator.

Defines:

- [`OUTPUT_FIELDS`](@ref): the names, in row order
- [`MonthlyAccumulator`](@ref): one field per row
- `accumulate!(acc, args...)`: adds one timestep, in a single `@turbo` loop
- `monthly_means(acc, nstep)`: the sums divided by `nstep`, as a record
- `reset!(acc)`: zeroes the sums
"""
macro output_fields(args, rows)
    names = Symbol[]
    arrays = Pair{Symbol,Expr}[]
    updates = Expr[]

    # `struct.field` becomes a local array indexed at the cell
    function at_cell(ex)
        ex isa Expr || return ex
        if ex.head === :.
            array = Symbol(ex.args[1], :_, ex.args[2].value)
            any(p -> first(p) === array, arrays) || push!(arrays, array => ex)
            return :($array[i, j])
        end
        return Expr(ex.head, map(at_cell, ex.args)...)
    end

    for row in rows.args
        row isa LineNumberNode && continue
        row isa Expr && row.head in (:+=, :-=) && row.args[1] isa Symbol ||
            error("@output_fields: expected `name += expr` or `name -= expr`, got `$row`")
        name = row.args[1]
        push!(names, name)
        push!(updates, Expr(row.head, :($(Symbol(:sum_, name))[i, j]), at_cell(row.args[2])))
    end

    return esc(quote
        const OUTPUT_FIELDS = $(Expr(:tuple, QuoteNode.(names)...))

        struct MonthlyAccumulator
            $((:($name::Matrix{Float32}) for name in names)...)
        end

        function accumulate!(acc::MonthlyAccumulator, $(args.args...))
            $((:($(Symbol(:sum_, name)) = acc.$name) for name in names)...)
            $((:($array = $field) for (array, field) in arrays)...)
            @turbo for j in 1:ydim
                for i in 1:xdim
                    $(updates...)
                end
            end
            return nothing
        end

        monthly_means(acc::MonthlyAccumulator, nstep) =
            $(Expr(:tuple, (Expr(:(=), name, :(acc.$name ./ nstep)) for name in names)...))

        function reset!(acc::MonthlyAccumulator)
            $((:(fill!(acc.$name, 0.0f0)) for name in names)...)
            return acc
        end
    end)
end

# `tend` is what `tendencies!` returns
@output_fields (surf::SurfaceState, tend, ws::ModelWorkspace) begin
    Ts     += surf.Ts
    Ta     += surf.Ta
    To     += surf.To
    q      += surf.q
    albedo += tend.albedo
    ice    += tend.ice_cover
    precip += ws.precip
    evap   += ws.evap
    qcrcl  += ws.qcrcl
    sw     += tend.SW
    lw     += tend.LW_surf
    qlat   += tend.Q_lat
    qsens  += tend.Q_sens
    olr    -= tend.LW_up + (1.0f0 - tend.em) * tend.LW_surf
    lwdown -= tend.LW_down
end

"""
    OUTPUT_FIELDS

The names of the monthly output fields, in the order of a
[`MonthlyRecord`](@ref). The record type and the accumulator are built from
the same list.
"""
OUTPUT_FIELDS

"""
    MonthlyAccumulator

The sums of the output fields over the current month, one `(xdim, ydim)` field
per name in [`OUTPUT_FIELDS`](@ref); `olr` and `lwdown` carry the signs of
[`MonthlyRecord`](@ref). [`output!`](@ref) adds to it every step, divides it
into a record at the end of the month and resets it with `reset!`.
"""
MonthlyAccumulator

MonthlyAccumulator() = MonthlyAccumulator(map(_ -> zeros(Float32, xdim, ydim), OUTPUT_FIELDS)...)

"""
    MonthlyRecord

One monthly-mean output record: a `NamedTuple` with fields `Ts`, `Ta`, `To`,
`q`, `albedo`, `ice`, `precip`, `evap`, `qcrcl`, `sw`, `lw`, `qlat`, `qsens`,
`olr`, `lwdown` ([`OUTPUT_FIELDS`](@ref)), each an `(xdim, ydim)`
`Matrix{Float32}`. Produced by [`output!`](@ref); `greb_model!`'s
`ctrl`/`scnr` results are `Vector{MonthlyRecord}`.

The fluxes are in W/m2. `sw`, `lw`, `qlat` and `qsens` are positive into the
surface; `lw` is the surface's own emission alone, so it is negative. `olr` is
the longwave leaving to space, positive upward; `lwdown` is the longwave the
air sends down, positive into the surface.
"""
const MonthlyRecord = NamedTuple{OUTPUT_FIELDS,NTuple{length(OUTPUT_FIELDS),Matrix{Float32}}}

"The climatologies a [`BoundaryAnomaly`](@ref) adds its anomalies to at scenario start."
const _BOUNDARY_FIELDS = (:Ts_clim, :u_clim, :v_clim, :omega_clim, :wind_speed_clim)

"""
    BoundaryAnomalyFields

The anomalies of a [`BoundaryAnomaly`](@ref) scenario: a `NamedTuple` with one
`(xdim, ydim, nstep_yr)` array per climatology the anomaly is added to
(`Ts_clim`, `u_clim`, `v_clim`, `omega_clim`, `wind_speed_clim`).
"""
const BoundaryAnomalyFields = NamedTuple{_BOUNDARY_FIELDS,NTuple{length(_BOUNDARY_FIELDS),Array{Float32,3}}}

"""
    ClimateFields

Loaded climatology, derived grid fields, flux corrections, the CO2 mask and
the insolation table: what `load_climatology` fills in and every physics
function reads. It is passed as an argument, never held as global state. One
loaded instance can serve several runs one after another: [`greb_model!`](@ref)
restores the input fields it changes and derives the others again at the
start of each run.

`ClimateFields()` builds an all-zero instance with `loaded = false`;
`load_climatology` sets `loaded = true`. [`greb_model!`](@ref) refuses unloaded
fields unless `allow_uninitialized=true`: an all-zero climatology raises no
error, it runs and returns NaN in every output field.

`boundary_anomaly` is `nothing` until [`load_boundary_anomaly!`](@ref) loads
the anomalies of a [`BoundaryAnomaly`](@ref) scenario; it then holds the
[`BoundaryAnomalyFields`](@ref) of one source at a time.
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
    Ts_clim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # surface temperature [K]
    u_clim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # zonal wind [m/s]
    v_clim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # meridional wind [m/s]
    q_clim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # atmospheric humidity [kg/kg]
    mld_clim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # mixed-layer depth [m]
    omega_clim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # vertical velocity [Pa/s]
    omega_std_clim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # omega std deviation [Pa/s]
    wind_speed_clim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # wind speed [m/s]

    # The winds split by sign (derive_fields!): pos + neg = the wind
    u_clim_pos::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # eastward part, 0 elsewhere
    u_clim_neg::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # westward part, 0 elsewhere
    v_clim_pos::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # northward part, 0 elsewhere
    v_clim_neg::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # southward part, 0 elsewhere

    To_clim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # deep ocean temperature [K]
    cloud_clim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # cloud cover fraction
    soil_wetness_clim::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # soil wetness [0-1]

    # Solar / radiation
    sw_solar::Matrix{Float32} = zeros(Float32, ydim, nstep_yr)  # 24hr mean solar radiation [W/m²] (ydim, nstep_yr)
    dTrad::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)  # Tatmos-radiation offset

    # Flux correction arrays (zeros unless loaded from file)
    Ts_flux_correction::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    q_flux_correction::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    To_flux_correction::Array{Float32,3} = zeros(Float32, xdim, ydim, nstep_yr)
    # File the three arrays were loaded from by `load_flux_corrections!`; "" if unknown.
    # A user who writes into the arrays sets it back to "".
    flux_source::String = ""

    # Regional CO2 mask (1.0 = full CO2, 0.5 = half CO2)
    co2_part::Matrix{Float32} = ones(Float32, xdim, ydim)

    # The anomalies of a BoundaryAnomaly scenario and where they were loaded
    # from (directory and source). Set by `load_boundary_anomaly!`.
    boundary_anomaly::Union{Nothing,BoundaryAnomalyFields} = nothing
    boundary_anomaly_source::Tuple{String,Symbol} = ("", :none)

    # false for a bare `ClimateFields()`; set by `load_climatology`. See the
    # docstring above and `greb_model!`'s `allow_uninitialized` keyword.
    loaded::Bool = false
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
radiation-temperature offset `dTrad` (from `Ts_clim`), the deep-ocean depth
`z_ocean` (from `mld_clim`), the pressure weights `wz_air`/`wz_vapor` and the
rain limit (from `z_topo`), the surface heat capacity `cap_surf` (from
`z_topo` and `mld_clim`) and the winds split by sign. Call it again after
changing an input map; [`init_model!`](@ref) calls it once per run.
"""
function derive_fields!(fields::ClimateFields, p::Processes)
    z_topo = fields.z_topo
    mld_clim = fields.mld_clim

    @. fields.dTrad = -0.16f0 * fields.Ts_clim - 5.0f0
    # Three times the deepest mixed layer of the year
    fields.z_ocean .= 3.0f0 .* dropdims(maximum(mld_clim; dims=3); dims=3)

    @. fields.wz_air = exp(-z_topo / z_air)
    @. fields.wz_vapor = exp(-z_topo / z_vapor)
    @. fields.rain_limit = -0.0015f0 / (fields.wz_vapor * r_qviwv * 86400.0f0)

    cap_surf = fields.cap_surf
    for j in 1:ydim, i in 1:xdim
        cap_surf[i, j] = is_land(z_topo[i, j]) || p.ocean === :none ? cap_land : cap_ocean * mld_clim[i, j, 1]
    end

    split_winds!(fields)
    return fields
end

# Upwind advection reads the eastward/westward and northward/southward parts
function split_winds!(fields::ClimateFields)
    _split_sign!(fields.u_clim_pos, fields.u_clim_neg, fields.u_clim)
    _split_sign!(fields.v_clim_pos, fields.v_clim_neg, fields.v_clim)
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

This is scratch space for the printed summary, not an output path - `Ts_annual_mean` is
averaged, printed and zeroed within a single `diagnostics!` call, so it never
holds a readable annual mean once the call returns. Model output is the
`Vector{MonthlyRecord}` that [`greb_model!`](@ref) returns.
"""
mutable struct ModelState
    sw_solar_forcing::Float32   # runtime solar multiplier used by SWradiation!
    Ts_annual_mean::Matrix{Float32}       # surface-temperature accumulator for the progress line
end

function ModelState()
    ModelState(1.0f0, zeros(Float32, xdim, ydim))
end
