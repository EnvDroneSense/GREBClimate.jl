# Which physics runs (Processes), the hydrology scheme (Hydrology) and the flux
# corrections (Corrections): three of the parts of a `Config`.

function _check_option(field, value, allowed)
    value in allowed || throw(ArgumentError("$field = :$value is not one of $(join((":$a" for a in allowed), ", "))"))
    return value
end

"""
    Processes(; kwargs...)

Which physical processes run. Every option says what the model does instead
when a process is off. The default is everything on, with the observed
climatologies.

| Keyword | Values (default first) | Meaning of the other values |
|:--------|:-----------------------|:----------------------------|
| `atmosphere` | `true`, `false` | No atmosphere: no longwave back-radiation, sensible or latent heat, transport |
| `clouds` | `:observed`, `:none`, `:uniform` | Cloud cover 0, or 0.7 everywhere |
| `humidity` | `:observed`, `:uniform` | Humidity climatology 0.0052 kg/kg everywhere |
| `hydrology` | `:full`, `:no_evap_rain`, `:none` | No evaporation and rain; or no water cycle at all (humidity 0, not updated) |
| `ocean` | `:full`, `:mixed_layer`, `:none` | A 50 m mixed layer without deep ocean; or land heat capacity everywhere |
| `topography` | `:observed`, `:flat` | Topography capped at 1 m |
| `control_co2` | `true`, `false` | Control CO₂ 0 ppm |
| `ice_albedo` | `true`, `false` | No ice-albedo feedback |
| `transport` | `true`, `false` | No atmospheric heat and moisture transport |
| `heat_diffusion`, `heat_advection`, `vapour_diffusion`, `vapour_advection`, `moisture_convergence` | `true`, `false` | The single transport term switched off |
"""
struct Processes
    atmosphere::Bool
    clouds::Symbol
    humidity::Symbol
    hydrology::Symbol
    ocean::Symbol
    topography::Symbol
    control_co2::Bool
    ice_albedo::Bool
    transport::Bool
    heat_diffusion::Bool
    heat_advection::Bool
    vapour_diffusion::Bool
    vapour_advection::Bool
    moisture_convergence::Bool
end

function Processes(; atmosphere=true, clouds=:observed, humidity=:observed, hydrology=:full,
    ocean=:full, topography=:observed, control_co2=true, ice_albedo=true, transport=true,
    heat_diffusion=true, heat_advection=true, vapour_diffusion=true, vapour_advection=true,
    moisture_convergence=true)
    return Processes(atmosphere,
        _check_option(:clouds, clouds, (:observed, :none, :uniform)),
        _check_option(:humidity, humidity, (:observed, :uniform)),
        _check_option(:hydrology, hydrology, (:full, :no_evap_rain, :none)),
        _check_option(:ocean, ocean, (:full, :mixed_layer, :none)),
        _check_option(:topography, topography, (:observed, :flat)),
        control_co2, ice_albedo, transport, heat_diffusion, heat_advection,
        vapour_diffusion, vapour_advection, moisture_convergence)
end

"""
    Hydrology(; rain=:fitted, rain_fit=:era, evaporation=:original)

The evaporation and rain scheme. The default is the fitted scheme of Stassen
et al. (2019); [`mscm_hydrology`](@ref) gives the original GREB scheme.

| Keyword | Values |
|:--------|:-------|
| `rain` | `:fitted` (Stassen et al. 2019), `:original` (proportional to humidity), `:rh` (+ relative humidity, with a rain-rate limit), `:omega` (+ vertical velocity), `:rh_omega` |
| `rain_fit` | `:era`, `:ncep`: which reanalysis the `:fitted` coefficients were fitted to; only for `rain = :fitted` |
| `evaporation` | `:original` (climatological wind plus a fixed gust term), `:skin` (skin temperature, land/ocean exchange coefficients), `:original_gust`, `:skin_gust` (the same with modified gust terms) |
"""
struct Hydrology
    rain::Symbol
    rain_fit::Symbol
    evaporation::Symbol
end

function Hydrology(; rain=:fitted, rain_fit=:era, evaporation=:original)
    _check_option(:rain, rain, (:fitted, :original, :rh, :omega, :rh_omega))
    _check_option(:rain_fit, rain_fit, (:era, :ncep))
    rain_fit === :ncep && rain !== :fitted &&
        throw(ArgumentError("rain_fit = :ncep only applies to rain = :fitted"))
    _check_option(:evaporation, evaporation, (:original, :skin, :original_gust, :skin_gust))
    return Hydrology(rain, rain_fit, evaporation)
end

"""
    mscm_hydrology()

The original GREB hydrology`Hydrology(rain = :original, evaporation = :original)`. 
Reproducing MSCM also needs `Processes(moisture_convergence = false)`.
"""
mscm_hydrology() = Hydrology(rain=:original, evaporation=:original)

"""
    Corrections

How the flux corrections that hold the control climate at the observed one
are obtained: [`SpinUp`](@ref)`(years)` computes them, [`Stored`](@ref)`()`
reads the dataset's file, [`NoCorrections`](@ref)`()` runs without.
"""
abstract type Corrections end

"""
    SpinUp(years)

Compute the flux corrections in a spin-up of `years` years before the control
run.
"""
struct SpinUp <: Corrections
    years::Int
    function SpinUp(years::Integer)
        years >= 0 || throw(ArgumentError("SpinUp years must be >= 0, got $years"))
        return new(years)
    end
end

"""
    Stored()

Use the flux corrections stored in the dataset instead of computing them.
"""
struct Stored <: Corrections end

"""
    NoCorrections()

Run without flux corrections; the control climate drifts from the observed one.
"""
struct NoCorrections <: Corrections end
