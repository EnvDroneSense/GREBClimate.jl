# What a scenario imposes: the CO2 path and where it applies, the sun, a surface
# forcing, the control CO2, the calendar start and the output form.

"""
    CO2Path

The scenario's CO₂ concentration over time: [`ConstantCO2`](@ref),
[`CO2Table`](@ref), [`CO2File`](@ref), [`A1BRamp`](@ref),
[`CO2SineWave`](@ref), [`CO2Step`](@ref), [`SeasonalCO2`](@ref).
"""
abstract type CO2Path end

"""
    ConstantCO2(ppm)

The same CO₂ concentration for the whole scenario.
"""
struct ConstantCO2 <: CO2Path
    ppm::Float32
end

"""
    CO2Table(key)

CO₂ per year from the dataset's IPCC table `key`: `:rcp26`, `:rcp45`, `:rcp6`,
`:ssp119`, `:ssp126`, `:ssp245`, `:ssp460`, `:ssp585` or `:hist` (1850-2017).
"""
struct CO2Table <: CO2Path
    key::Symbol
end

"""
    CO2File(path)

CO₂ per year from a text file of `year ppm` lines (see
[`load_custom_co2_scenario`](@ref)).
"""
struct CO2File <: CO2Path
    path::String
end

"""
    A1BRamp()

The SRES A1B path: 310 ppm in 1950, 370 in 2000, 520 in 2050, 700 in 2100,
linear in between; 340 ppm after 2100.
"""
struct A1BRamp <: CO2Path end

"""
    CO2SineWave()

340 to 680 ppm in a 30-year cosine cycle: `510 + 170 cos(2π (year - 13) / 30)`.
"""
struct CO2SineWave <: CO2Path end

"""
    CO2Step(before, after, year)

`before` ppm until `year`, `after` ppm from `year` on.
"""
struct CO2Step <: CO2Path
    before::Float32
    after::Float32
    year::Int
end

"""
    SeasonalCO2(inside, outside, season)

`inside` ppm during `season` (`:boreal_winter`, steps 1-181 and 547-730 of the
year, or `:boreal_summer`, the rest), `outside` ppm otherwise.
"""
struct SeasonalCO2 <: CO2Path
    inside::Float32
    outside::Float32
    season::Symbol
    function SeasonalCO2(inside, outside, season::Symbol)
        _check_option(:season, season, (:boreal_winter, :boreal_summer))
        return new(inside, outside, season)
    end
end

"""
    CO2Mask

Where the scenario CO₂ applies; elsewhere the cell gets half of it (the
regional 2×CO₂ experiments): [`UniformMask`](@ref), [`LatitudeMask`](@ref),
[`SurfaceMask`](@ref).
"""
abstract type CO2Mask end

"""
    UniformMask()

The scenario CO₂ applies everywhere.
"""
struct UniformMask <: CO2Mask end

"""
    LatitudeMask(band)

The scenario CO₂ applies in `band`: `:nh`, `:sh`, `:tropics` or `:extratropics`.
"""
struct LatitudeMask <: CO2Mask
    band::Symbol
    function LatitudeMask(band::Symbol)
        _check_option(:band, band, (:nh, :sh, :tropics, :extratropics))
        return new(band)
    end
end

"""
    SurfaceMask(surface)

The scenario CO₂ applies over `:ocean` or over `:land_ice`, as the control
run's annual-mean ice cover defines them.
"""
struct SurfaceMask <: CO2Mask
    surface::Symbol
    function SurfaceMask(surface::Symbol)
        _check_option(:surface, surface, (:ocean, :land_ice))
        return new(surface)
    end
end

"""
    Solar

The incoming sunlight: [`ModernSolar`](@ref), [`SolarConstant`](@ref),
[`SolarCycle`](@ref), [`SolarTable`](@ref), [`EarthSunDistance`](@ref).
"""
abstract type Solar end

"""
    ModernSolar()

Today's insolation from the dataset.
"""
struct ModernSolar <: Solar end

"""
    SolarConstant(dW)

The solar constant raised from 1365 by `dW` W/m².
"""
struct SolarConstant <: Solar
    dW::Float32
end

"""
    SolarCycle(amplitude, period)

The solar constant varying as `1365 + amplitude sin(2π year / period)` W/m².
"""
struct SolarCycle <: Solar
    amplitude::Float32
    period::Float32
end

"""
    SolarTable(kind, index)

Insolation from the dataset's solar scenario tables: `kind` is `:paleo` (231
kyr ago), `:obliquity` or `:eccentricity`; `index` selects the table row of
the last two. Obliquity row `k` is `-25 + k/2` degrees (`k` = 0, 5, ..., 230),
eccentricity row `k` is `-0.30 + 0.01k` (`k` = 0, ..., 60; negative puts
perihelion in July). The default is the row nearest today: 95 (22.5 degrees)
and 32 (0.02).
"""
struct SolarTable <: Solar
    kind::Symbol
    index::Int
    function SolarTable(kind::Symbol, index::Integer=get(_MODERN_ROW, kind, 0))
        _check_option(:kind, kind, (:paleo, :obliquity, :eccentricity))
        return new(kind, index)
    end
end

const _MODERN_ROW = Dict(:obliquity => 95, :eccentricity => 32)

"""
    EarthSunDistance(pct)

The Earth-Sun distance changed by `pct` percent; insolation scales with its
inverse square.
"""
struct EarthSunDistance <: Solar
    pct::Float32
end

"""
    SurfaceForcing

A forcing imposed on the surface during the scenario:
[`NoSurfaceForcing`](@ref), [`BoundaryAnomaly`](@ref), [`SSTOffset`](@ref).
"""
abstract type SurfaceForcing end

"""
    NoSurfaceForcing()

No imposed surface forcing.
"""
struct NoSurfaceForcing <: SurfaceForcing end

"""
    BoundaryAnomaly(source)

Anomalies of surface temperature, winds and vertical velocity added to the
climatology for the scenario, with surface temperature held at it: `source`
is `:cmip5_rcp85`, `:elnino` or `:lanina`.
"""
struct BoundaryAnomaly <: SurfaceForcing
    source::Symbol
    function BoundaryAnomaly(source::Symbol)
        _check_option(:source, source, (:cmip5_rcp85, :elnino, :lanina))
        return new(source)
    end
end

"""
    SSTOffset(K)

Ocean surface temperature held at the climatology plus `K`, with CO₂ at its
control value.
"""
struct SSTOffset <: SurfaceForcing
    K::Float32
end

"""
    Scenario(; co2=ConstantCO2(340), co2_mask=UniformMask(), solar=ModernSolar(),
               surface=NoSurfaceForcing(), control_co2=340, start_year=1950, output=:anomaly)

What an experiment imposes. `control_co2` is the CO₂ of the spin-up and the
control run; `co2` the scenario's. `start_year` is the calendar year the
scenario starts in; `output` is `:anomaly` (each month minus the control's
final year) or `:absolute`.
"""
struct Scenario
    co2::CO2Path
    co2_mask::CO2Mask
    solar::Solar
    surface::SurfaceForcing
    control_co2::Float32
    start_year::Int
    output::Symbol
end

function Scenario(; co2=ConstantCO2(340), co2_mask=UniformMask(), solar=ModernSolar(),
    surface=NoSurfaceForcing(), control_co2=340, start_year=1950, output=:anomaly)
    _check_option(:output, output, (:anomaly, :absolute))
    return Scenario(co2, co2_mask, solar, surface, control_co2, start_year, output)
end
