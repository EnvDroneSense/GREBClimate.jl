# A full model configuration, and the named experiments as presets of it.

"""
    Config(; scenario=Scenario(), processes=Processes(), hydrology=Hydrology(),
             corrections=SpinUp(3), modules=())

Everything a run needs besides its length ([`RunSpec`](@ref)) and the data:
what the experiment imposes ([`Scenario`](@ref)), which physics runs
([`Processes`](@ref)), the hydrology scheme ([`Hydrology`](@ref)), how the flux
corrections are obtained ([`Corrections`](@ref)) and optional add-on modules
(none yet). The default is the full model under constant 340 ppm. Build named
experiments with [`preset`](@ref).
"""
Base.@kwdef struct Config
    scenario::Scenario = Scenario()
    processes::Processes = Processes()
    hydrology::Hydrology = Hydrology()
    corrections::Corrections = SpinUp(3)
    modules::Tuple = ()
end

# The scenario of every named experiment. IPCC-table and A1B runs start from a
# 280 ppm control; orbital runs start in year 1 and are returned absolute.
const _TABLE = (control_co2=280,)
const _ORBITAL = (start_year=1, output=:absolute)
const _PRESET_SCENARIOS = Dict{Symbol,Scenario}(
    :full_model => Scenario(),
    :co2_double => Scenario(co2=ConstantCO2(680)),
    :co2_quadruple => Scenario(co2=ConstantCO2(1360)),
    :co2_10x => Scenario(co2=ConstantCO2(3400)),
    :co2_half => Scenario(co2=ConstantCO2(170)),
    :co2_zero => Scenario(co2=ConstantCO2(0)),
    :co2_sine_wave => Scenario(co2=CO2SineWave()),
    :co2_abrupt_reverse => Scenario(co2=CO2Step(680, 340, 1980)),
    :a1b => Scenario(; co2=A1BRamp(), _TABLE...),
    :solar_plus27 => Scenario(solar=SolarConstant(27)),
    :solar_cycle_11yr => Scenario(solar=SolarCycle(1, 11)),
    :paleo_231kyr => Scenario(co2=ConstantCO2(200), solar=SolarTable(:paleo)),
    :paleo_solar_modern_co2 => Scenario(solar=SolarTable(:paleo)),
    :modern_solar_paleo_co2 => Scenario(co2=ConstantCO2(200)),
    :obliquity => Scenario(; solar=SolarTable(:obliquity), _ORBITAL...),
    :eccentricity => Scenario(; solar=SolarTable(:eccentricity), _ORBITAL...),
    :earth_sun_distance => Scenario(; solar=EarthSunDistance(0), _ORBITAL...),
    :elnino => Scenario(surface=BoundaryAnomaly(:elnino)),
    :lanina => Scenario(surface=BoundaryAnomaly(:lanina)),
    :rcp85_boundary => Scenario(surface=BoundaryAnomaly(:cmip5_rcp85)),
    :rcp26 => Scenario(; co2=CO2Table(:rcp26), _TABLE...),
    :rcp45 => Scenario(; co2=CO2Table(:rcp45), _TABLE...),
    :rcp60 => Scenario(; co2=CO2Table(:rcp6), _TABLE...),
    :ssp119 => Scenario(; co2=CO2Table(:ssp119), _TABLE...),
    :ssp126 => Scenario(; co2=CO2Table(:ssp126), _TABLE...),
    :ssp245 => Scenario(; co2=CO2Table(:ssp245), _TABLE...),
    :ssp460 => Scenario(; co2=CO2Table(:ssp460), _TABLE...),
    :ssp585 => Scenario(; co2=CO2Table(:ssp585), _TABLE...),
    :historical_co2 => Scenario(; co2=CO2Table(:hist), start_year=1850, _TABLE...),
    :custom_co2 => Scenario(; co2=CO2File(""), _TABLE...),
    :sst_plus1 => Scenario(surface=SSTOffset(1)),
    :regional_co2_nh => Scenario(co2=ConstantCO2(680), co2_mask=LatitudeMask(:nh)),
    :regional_co2_sh => Scenario(co2=ConstantCO2(680), co2_mask=LatitudeMask(:sh)),
    :regional_co2_tropics => Scenario(co2=ConstantCO2(680), co2_mask=LatitudeMask(:tropics)),
    :regional_co2_extratropics => Scenario(co2=ConstantCO2(680), co2_mask=LatitudeMask(:extratropics)),
    :regional_co2_ocean => Scenario(co2=ConstantCO2(680), co2_mask=SurfaceMask(:ocean)),
    :regional_co2_land_ice => Scenario(co2=ConstantCO2(680), co2_mask=SurfaceMask(:land_ice)),
    :regional_co2_winter => Scenario(co2=SeasonalCO2(680, 340, :boreal_winter)),
    :regional_co2_summer => Scenario(co2=SeasonalCO2(680, 340, :boreal_summer)),
    # Deconstruction experiments: switch processes off with `processes`
    :decon_mean_climate => Scenario(),
    :decon_2xco2 => Scenario(co2=ConstantCO2(680)),
)

"""
    preset_names()

The names [`preset`](@ref) accepts, sorted.
"""
preset_names() = sort!(collect(keys(_PRESET_SCENARIOS)))

"""
    preset(name; processes=Processes(), hydrology=Hydrology(), corrections=nothing,
           index=0, pct=0, path="") -> Config

The [`Config`](@ref) of a named experiment ([`preset_names`](@ref) lists them).
`processes` and `hydrology` change the physics, typically for
`:decon_mean_climate` (mean-climate deconstruction; run it with
`RunSpec(scnr = 0)`) and `:decon_2xco2` (2×CO₂-response deconstruction).
`index` selects the table row of `:obliquity`/`:eccentricity`, `pct` the
distance change of `:earth_sun_distance`, `path` the CO₂ file of
`:custom_co2`. `corrections` defaults to a 3-year [`SpinUp`](@ref), or to the
[`Stored`](@ref) corrections when topography is flat.

```jldoctest
julia> preset(:co2_double).scenario.co2
ConstantCO2(680.0f0)
```
"""
function preset(name::Symbol; processes::Processes=Processes(), hydrology::Hydrology=Hydrology(),
    corrections::Union{Corrections,Nothing}=nothing, index::Integer=0, pct::Real=0, path::AbstractString="")
    haskey(_PRESET_SCENARIOS, name) ||
        throw(ArgumentError("unknown preset :$name; known: $(join((":$p" for p in preset_names()), ", "))"))
    s = _PRESET_SCENARIOS[name]
    s.solar isa SolarTable && s.solar.kind !== :paleo && (s = _with(s; solar=SolarTable(s.solar.kind, index)))
    s.solar isa EarthSunDistance && (s = _with(s; solar=EarthSunDistance(pct)))
    s.co2 isa CO2File && (s = _with(s; co2=CO2File(path)))
    corrections = something(corrections, processes.topography === :flat ? Stored() : SpinUp(3))
    return Config(; scenario=s, processes, hydrology, corrections)
end

_with(s::Scenario; kw...) = Scenario(; (f => getfield(s, f) for f in fieldnames(Scenario))..., kw...)
