# The net heat flux into the surface and into the atmosphere (W/m2), defined
# once for the flux-correction spin-up and the run. Macros, not functions: a
# function call changes how `@turbo` contracts `em * LW_surf` and so changes
# `Ta` in the last digit. They read the flux arrays by name from the caller.

macro surface_flux(i, j)
    return esc(:(SW[$i, $j] + LW_surf[$i, $j] - LW_down[$i, $j] + Q_lat[$i, $j] + Q_sens[$i, $j]))
end

macro atmosphere_flux(i, j)
    return esc(:(LW_up[$i, $j] + LW_down[$i, $j] - em[$i, $j] * LW_surf[$i, $j] + Q_lat_air[$i, $j] - Q_sens[$i, $j]))
end

# The same sums at one cell of a `tendencies!` result, for tests and tools.
# Outside `@turbo` they can differ from the loops' values in the last digit.
function surface_flux(tend, i, j)
    SW, LW_surf, LW_down, Q_lat, Q_sens = tend.SW, tend.LW_surf, tend.LW_down, tend.Q_lat, tend.Q_sens
    return @surface_flux(i, j)
end

function atmosphere_flux(tend, i, j)
    LW_up, LW_down, em, LW_surf, Q_lat_air, Q_sens = tend.LW_up, tend.LW_down, tend.em, tend.LW_surf, tend.Q_lat_air, tend.Q_sens
    return @atmosphere_flux(i, j)
end

# What an observer of `greb_model!` is handed at each call. The arrays are the
# model's own, not copies.
_step_view(t::ModelTime, CO2, Ts, Ta, To, q, tend, fields, r) =
    (; time=t, phase=t.phase, step=t.step, year=year(t), step_of_year=step_of_year(t),
        CO2, Ts, Ta, To, q, tend, fields, config=r)

"""
    BudgetCheck()

An observer for [`greb_model!`](@ref)'s `observer` keyword that checks the
model's bookkeeping. Each step it works out what the surface, air and ocean
temperatures and the humidity should become from the step's flows (the
limiters included) and compares that with what the model stored.

After the run it holds, over all steps and cells:

| Field | Meaning |
|:------|:--------|
| `steps` | steps checked |
| `surface`, `atmosphere`, `ocean` | largest difference between the expected and the stored `Ts`, `Ta`, `To` [K] |
| `water` | the same for `q` [kg/kg] |
| `floor_Ts`, `floor_Ta` | cells held at the temperature floor `min_T_K` |
| `humidity_low`, `humidity_high` | cells where the humidity change was limited |
| `rain_limit` | cells where rain was set to the rain limit (`rain = :rh` only) |

Not exported: create it with `GREBClimate.BudgetCheck()`.
"""
Base.@kwdef mutable struct BudgetCheck
    Ts::Matrix{Float32} = zeros(Float32, xdim, ydim)
    Ta::Matrix{Float32} = zeros(Float32, xdim, ydim)
    To::Matrix{Float32} = zeros(Float32, xdim, ydim)
    q::Matrix{Float32} = zeros(Float32, xdim, ydim)
    cap_surf::Matrix{Float32} = zeros(Float32, xdim, ydim)
    steps::Int = 0
    surface::Float32 = 0.0f0
    atmosphere::Float32 = 0.0f0
    ocean::Float32 = 0.0f0
    water::Float32 = 0.0f0
    floor_Ts::Int = 0
    floor_Ta::Int = 0
    humidity_low::Int = 0
    humidity_high::Int = 0
    rain_limit::Int = 0
end

function (b::BudgetCheck)(point::Symbol, view)
    if point === :after_tendencies
        # The update overwrites the state, and seaice! the heat capacity
        b.Ts .= view.Ts; b.Ta .= view.Ta; b.To .= view.To; b.q .= view.q
        b.cap_surf .= view.fields.cap_surf
    elseif point === :after_step
        _check_step!(b, view)
    end
    return nothing
end

_floor(T) = ifelse(T < min_T_K, min_T_K, T)

function _check_step!(b::BudgetCheck, view)
    tend, fields, slice = view.tend, view.fields, view.step_of_year
    config = view.config
    hydro_on = config.config.processes.hydrology !== :none
    rain_limited = config.config.processes.hydrology === :full && config.config.processes.atmosphere &&
                   config.hydrology.rain === :rh
    for j in 1:ydim, i in 1:xdim
        Ts = b.Ts[i, j] + tend.dT_ocean[i, j] +
             Δt * (surface_flux(tend, i, j) + fields.Ts_flux_correction[i, j, slice]) / b.cap_surf[i, j]
        Ta = b.Ta[i, j] + tend.dTa_crcl[i, j] + Δt * atmosphere_flux(tend, i, j) / cap_air
        To = b.To[i, j] + tend.dTo[i, j] + fields.To_flux_correction[i, j, slice]
        dq = Δt * (tend.dq_eva[i, j] + tend.dq_rain[i, j]) + tend.dq_crcl[i, j] + fields.q_flux_correction[i, j, slice]
        low = dq <= -b.q[i, j]
        low && (dq = -min_humidity_change * b.q[i, j])
        high = dq > max_humidity_change
        high && (dq = max_humidity_change)

        b.floor_Ts += Ts < min_T_K
        b.floor_Ta += Ta < min_T_K
        b.humidity_low += hydro_on & low
        b.humidity_high += hydro_on & high
        b.rain_limit += rain_limited & (tend.dq_rain[i, j] == fields.rain_limit[i, j])

        b.surface = max(b.surface, abs(view.Ts[i, j] - _floor(Ts)))
        b.atmosphere = max(b.atmosphere, abs(view.Ta[i, j] - _floor(Ta)))
        b.ocean = max(b.ocean, abs(view.To[i, j] - To))
        b.water = max(b.water, abs(view.q[i, j] - (b.q[i, j] + hydro_on * dq)))
    end
    b.steps += 1
    return b
end

"""
    RangeCheck(; Ts=(100, 360), Ta=(100, 360), To=(250, 330), q=(0, 0.08))

An observer for [`greb_model!`](@ref)'s `observer` keyword that notices a run
which leaves the physical range without going non-finite: a cell can run away
to 1e24 K and still look like a result. Each keyword is the allowed
`(lowest, highest)` value of that field, in K and kg/kg; the defaults are well
outside what any preset reaches (`Ts` 155 to 325 K, `Ta` 162 to 329 K, `To` 266
to 313 K, `q` up to 0.039 kg/kg).

After the run it holds:

| Field | Meaning |
|:------|:--------|
| `steps` | steps checked |
| `seen` | lowest and highest value of each field over the run |
| `first` | `nothing`, or the first step outside the range: `(phase, step, year, field, value)` |

`GREBClimate.in_range(check)` is `true` when no step left the range. A NaN
counts as outside. Not exported: create it with `GREBClimate.RangeCheck()`.
"""
mutable struct RangeCheck
    limits::NamedTuple{(:Ts, :Ta, :To, :q),NTuple{4,Tuple{Float32,Float32}}}
    seen::NamedTuple{(:Ts, :Ta, :To, :q),NTuple{4,Tuple{Float32,Float32}}}
    steps::Int
    first::Union{Nothing,NamedTuple{(:phase, :step, :year, :field, :value),Tuple{Phase,Int,Int,Symbol,Float32}}}
end

function RangeCheck(; Ts=(100, 360), Ta=(100, 360), To=(250, 330), q=(0, 0.08))
    nothing_seen = (Inf32, -Inf32)
    return RangeCheck(map(l -> Float32.(l), (; Ts, Ta, To, q)), map(_ -> nothing_seen, (; Ts, Ta, To, q)), 0, nothing)
end

in_range(c::RangeCheck) = c.first === nothing

# Lowest and highest value of `A` and the number of cells outside
# `lower..upper`. A NaN fails both comparisons, so it counts as outside; the
# extremes skip it.
function _range(A::AbstractMatrix{Float32}, lower::Float32, upper::Float32)
    lo, hi, outside = Inf32, -Inf32, 0
    @turbo for i in eachindex(A)
        x = A[i]
        lo = min(lo, x)
        hi = max(hi, x)
        outside += !((lower <= x) & (x <= upper))
    end
    return lo, hi, outside
end

function (c::RangeCheck)(point::Symbol, view)
    point === :after_step || return nothing
    seen = map(keys(c.limits)) do field
        old = c.seen[field]
        lower, upper = c.limits[field]
        lo, hi, outside = _range(getproperty(view, field), lower, upper)
        if c.first === nothing && outside > 0
            value = lo < lower ? lo : hi > upper ? hi : NaN32
            c.first = (phase=view.phase, step=Int(view.step), year=Int(view.year), field=field, value=value)
        end
        (min(old[1], lo), max(old[2], hi))
    end
    c.seen = NamedTuple{keys(c.limits)}(seen)
    c.steps += 1
    return nothing
end
