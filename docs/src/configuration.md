# Configuration

```@meta
DocTestSetup = :(using GREBClimate)
```

A run is described by a [`Config`](@ref) and its length by a
[`RunSpec`](@ref). Build a `Config` from a named experiment with
[`preset`](@ref), or from its parts:

| Part | Type | Decides |
|:-----|:-----|:--------|
| `scenario` | [`Scenario`](@ref) | What the experiment imposes: CO₂ path and where it applies, sunlight, surface forcing, control CO₂, start year, output form |
| `processes` | [`Processes`](@ref) | Which physical processes run, each with a stated meaning when off |
| `hydrology` | [`Hydrology`](@ref) | The evaporation and rain scheme |
| `corrections` | [`Corrections`](@ref) | How the flux corrections are obtained: [`SpinUp`](@ref)`(years)`, [`Stored`](@ref)`()`, [`NoCorrections`](@ref)`()` |

```julia
config = preset(:co2_double)                                   # a named experiment
config = preset(:co2_double; processes = (clouds = :uniform,))  # with one process changed
config = Config(scenario = Scenario(co2 = ConstantCO2(500)))    # from parts
```

Every part validates its options when it is built, so a misspelt option fails
at once instead of being ignored. The option tables below are the types' own
documentation.

## Presets

[`preset_names`](@ref)`()` lists them. The table is generated from the
presets; an empty cell means the default of [`Scenario`](@ref) (340 ppm
everywhere, modern sunlight, no surface forcing, control at 340 ppm, start in
1950, output as an anomaly).

```@eval
using GREBClimate, Markdown
d = Scenario()
cell(x, default) = x == default ? "" : "`" * replace(repr(x), "GREBClimate." => "", r"(\d)f0" => s"\1") * "`"
physics(c) = join(filter(!isempty, [
    c.processes == Processes() ? "" : "processes changed",
    c.hydrology == Hydrology() ? "" : (c.hydrology == mscm_hydrology() ? "MSCM hydrology" : "hydrology changed"),
    c.corrections == SpinUp(3) ? "" : "`" * repr(c.corrections) * "`"]), ", ")
rows = map(preset_names()) do p
    c = preset(p)
    s = c.scenario
    "| `:$p` | " * join([cell(s.co2, d.co2), cell(s.co2_mask, d.co2_mask), cell(s.solar, d.solar),
                          cell(s.surface, d.surface), cell(s.control_co2, d.control_co2),
                          cell(s.start_year, d.start_year), cell(s.output, d.output), physics(c)], " | ") * " |"
end
Markdown.parse("| Preset | CO₂ | Where | Sunlight | Surface | Control CO₂ | Start | Output | Physics |\n" *
               "|:--|:--|:--|:--|:--|:--|:--|:--|:--|\n" * join(rows, "\n"))
```

Three presets take a parameter: `preset(:obliquity; index = ...)` and
`preset(:eccentricity; index = ...)` pick the table row (default: the row
nearest today), `preset(:earth_sun_distance; pct = ...)` the distance change,
`preset(:custom_co2; path = ...)` the CO₂ file. `:decon_mean_climate` is the
mean-climate deconstruction: switch processes off with
`processes = (...)` and run it with `RunSpec(scnr = 0)`. On its stored
corrections a switched-off process changes the climate; a spin-up would
recompute the corrections for each configuration and pull every one back to
the observed climate. `:decon_2xco2` is the 2×CO₂-response deconstruction.

## Writing a scenario from parts

A [`Scenario`](@ref) combines one part of each kind; any combination runs.

```jldoctest scenario
julia> s = Scenario(co2 = CO2Step(340, 680, 2000), solar = SolarCycle(2, 11), control_co2 = 340);

julia> r = resolve(Config(scenario = s));

julia> forcing(1, 1999, r).CO2, forcing(1, 2000, r).CO2
(340.0f0, 680.0f0)
```

| Kind | Options |
|:-----|:--------|
| CO₂ path ([`CO2Path`](@ref)) | [`ConstantCO2`](@ref), [`CO2Table`](@ref), [`CO2File`](@ref), [`A1BRamp`](@ref), [`CO2SineWave`](@ref), [`CO2Step`](@ref), [`SeasonalCO2`](@ref) |
| Where it applies ([`CO2Mask`](@ref)) | [`UniformMask`](@ref), [`LatitudeMask`](@ref), [`SurfaceMask`](@ref) |
| Sunlight ([`Solar`](@ref)) | [`ModernSolar`](@ref), [`SolarConstant`](@ref), [`SolarCycle`](@ref), [`SolarTable`](@ref), [`EarthSunDistance`](@ref) |
| Surface ([`SurfaceForcing`](@ref)) | [`NoSurfaceForcing`](@ref), [`BoundaryAnomaly`](@ref), [`SSTOffset`](@ref) |

A scenario that reads a table ([`CO2Table`](@ref), [`CO2File`](@ref),
[`SolarTable`](@ref)) loads it when the configuration is resolved:
[`greb_model!`](@ref) does this itself, from its `jld2_dir`.

## Configuration types

```@autodocs
Modules = [GREBClimate]
Pages = ["config/presets.jl", "config/processes.jl", "config/scenario.jl", "config/resolve.jl"]
```

## References

1. Dommenget, D., and Flöter, J. (2011). Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model. *Climate Dynamics*, 37: 2143. [doi:10.1007/s00382-011-1026-0](https://doi.org/10.1007/s00382-011-1026-0)
2. Stassen, C., Dommenget, D., and Loveday, N. (2019). A hydrological cycle model for the Globally Resolved Energy Balance (GREB) model v1.0. *Geoscientific Model Development*, 12, 425-440. [doi:10.5194/gmd-12-425-2019](https://doi.org/10.5194/gmd-12-425-2019)

See the repository [README](https://github.com/EnvDroneSense/GREBClimate.jl#references)'s References section for the original GREB model homepage link.
