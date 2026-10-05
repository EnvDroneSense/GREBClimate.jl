# A Config with what its names refer to: the rain coefficients of its
# hydrology scheme and the CO2 and solar tables its scenario reads.

"""
    ResolvedHydrology

A [`Hydrology`](@ref) scheme with its rain-regression coefficients `c_q`,
`c_rq`, `c_omega`, `c_omega_std`; what [`hydro!`](@ref) reads. Build it with
[`resolve`](@ref)`(hydrology)`.
"""
struct ResolvedHydrology
    rain::Symbol
    evaporation::Symbol
    c_q::Float32
    c_rq::Float32
    c_omega::Float32
    c_omega_std::Float32
end

# (c_q, c_rq, c_omega, c_omega_std) per (rain, rain_fit); only :fitted has an
# :ncep fit, the constructor of Hydrology rejects it for the others
const _RAIN_COEFFICIENTS = Dict(
    (:original, :era) => (1.0f0, 0.0f0, 0.0f0, 0.0f0),
    (:rh, :era) => (-1.391649f0, 3.018774f0, 0.0f0, 0.0f0),
    (:omega, :era) => (0.862162f0, 0.0f0, -29.02096f0, 0.0f0),
    (:rh_omega, :era) => (-0.2685845f0, 1.4591853f0, -26.9858807f0, 0.0f0),
    (:fitted, :era) => (-1.88f0, 2.25f0, -17.69f0, 59.07f0),
    (:fitted, :ncep) => (-1.27f0, 1.99f0, -16.54f0, 21.15f0),
)

"""
    ResolvedConfig

A [`Config`](@ref) with the data its names refer to: `config`, `hydrology`
(a [`ResolvedHydrology`](@ref)), `co2_table` (year => ppm, empty unless the
scenario's CO2 comes from a [`CO2Table`](@ref) or [`CO2File`](@ref)) and
`solar_table` (the [`SolarTable`](@ref) insolation, or `nothing`). Build it
with [`resolve`](@ref); [`greb_model!`](@ref) does so itself.
"""
struct ResolvedConfig
    config::Config
    hydrology::ResolvedHydrology
    co2_table::Dict{Int,Float32}
    solar_table::Union{Nothing,Matrix{Float32}}
end

"""
    resolve(config::Config; jld2_dir="") -> ResolvedConfig
    resolve(hydrology::Hydrology) -> ResolvedHydrology

Replace the names in `config` with what they refer to: the rain coefficients
of its hydrology scheme, and the CO2 or solar table its scenario reads from
the dataset in `jld2_dir` (or from the user's CO2 file).

```jldoctest
julia> resolve(mscm_hydrology()).c_q
1.0f0
```
"""
function resolve(config::Config; jld2_dir::AbstractString="")
    s = config.scenario
    return ResolvedConfig(config, resolve(config.hydrology), _co2_table(s.co2, jld2_dir),
                          _solar_table(s.solar, jld2_dir))
end

function resolve(h::Hydrology)
    return ResolvedHydrology(h.rain, h.evaporation, _RAIN_COEFFICIENTS[(h.rain, h.rain_fit)]...)
end

_co2_table(::CO2Path, jld2_dir) = Dict{Int,Float32}()
_co2_table(c::CO2Table, jld2_dir) = load_co2_scenario(String(jld2_dir), c.key)
function _co2_table(c::CO2File, jld2_dir)
    isempty(c.path) && throw(ArgumentError("CO2File needs the path of a CO2 file"))
    return load_co2_custom(c.path)
end

_solar_table(::Solar, jld2_dir) = nothing
_solar_table(s::SolarTable, jld2_dir) = Matrix{Float32}(load_solar_forcing(String(jld2_dir), s.kind, s.index))
