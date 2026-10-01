# The net heat flux into the surface and into the atmosphere (W/m2): the sums
# both the flux-correction spin-up and the run integrate, defined once so a new
# term enters both.
#
# They are macros, not functions: a macro pastes the sum into each `@turbo`
# loop, where a function call changes how the loop contracts `em * LW_surf`
# into a fused multiply-add and so changes `Ta` in the last digit. The macros
# read the flux arrays by name from the calling scope: `SW`, `LW_surf`,
# `LW_down`, `Q_lat`, `Q_sens` for the surface; `LW_up`, `LW_down`, `em`,
# `LW_surf`, `Q_lat_air`, `Q_sens` for the atmosphere.

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
_step_view(phase, it, year, ityr, CO2, Ts, Ta, To, q, tend, fields, r) =
    (; phase, it, year, ityr, CO2, Ts, Ta, To, q, tend, fields, config=r)

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

Experimental and not exported, like the `observer` keyword itself.
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
    tend, fields, ityr = view.tend, view.fields, view.ityr
    config = view.config
    hydro_on = config.config.processes.hydrology !== :none
    rain_limited = config.config.processes.hydrology === :full && config.config.processes.atmosphere &&
                   config.hydrology.rain === :rh
    for j in 1:ydim, i in 1:xdim
        Ts = b.Ts[i, j] + tend.dT_ocean[i, j] +
             Δt * (surface_flux(tend, i, j) + fields.TF_correct[i, j, ityr]) / b.cap_surf[i, j]
        Ta = b.Ta[i, j] + tend.dTa_crcl[i, j] + Δt * atmosphere_flux(tend, i, j) / cap_air
        To = b.To[i, j] + tend.dTo[i, j] + fields.ToF_correct[i, j, ityr]
        dq = Δt * (tend.dq_eva[i, j] + tend.dq_rain[i, j]) + tend.dq_crcl[i, j] + fields.qF_correct[i, j, ityr]
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
