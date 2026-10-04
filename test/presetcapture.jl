# What each experiment preset is: the forcing it imposes and the physics it
# runs, captured without the dataset. Used by test_presets.jl and by
# tools/validation/preset_reference.jl, which writes the reference in
# test/data/preset_reference.jl. Needs testutils.jl.

const PRESETS = preset_names()
const SAMPLE_STEPS = (1, 200, 400, 600)   # step of the year: winter, summer, summer, winter
const NYEARS = 150

# Run-length encoding: [(value, count), ...]
function rle(v)
    runs = Tuple{eltype(v),Int}[]
    for x in v
        if !isempty(runs) && isequal(runs[end][1], x)
            runs[end] = (x, runs[end][2] + 1)
        else
            push!(runs, (x, 1))
        end
    end
    return runs
end
unrle(runs) = reduce(vcat, [fill(x, n) for (x, n) in runs])

_table_name(::CO2Path) = :none
_table_name(c::CO2Table) = c.key
_table_name(::CO2File) = :custom
_solar_name(::Solar) = :none
_solar_name(s::SolarTable) = s.kind
_boundary_name(::SurfaceForcing) = :none
_boundary_name(b::BoundaryAnomaly) = b.source === :cmip5_rcp85 ? :rcp85 : b.source
_boundary_name(::SSTOffset) = :sst_plus1

# A config part as field = value pairs, so the reference still parses when a
# part gains a field
_parts(x) = NamedTuple{fieldnames(typeof(x))}(Tuple(getfield(x, f) for f in fieldnames(typeof(x))))

function capture_preset(p::Symbol)
    config = preset(p)
    s = config.scenario
    co2_table = _table_name(s.co2)
    # stand-in for the table the run loads
    table = co2_table === :none ? Dict{Int,Float32}() :
            Dict(y => 280.0f0 + 0.5f0 * (y - 1850) for y in s.start_year:(s.start_year + NYEARS))
    r = ResolvedConfig(config, resolve(config.hydrology), table, nothing)
    fields = synthetic_fields()
    ini = quiet() do
        init_model!(r, fields)
    end
    # The two mask steps of the scenario start
    GREBClimate.apply_co2_mask!(s.co2_mask, fields)
    static_mask = copy(fields.co2_part)
    ice = zeros(Float32, X, Y, 12)
    ice[:, abs.(GREBClimate.lat_grid) .> 60, :] .= 1.0f0
    GREBClimate.apply_dynamic_co2_mask!(s.co2_mask, fields, ice)

    co2, solar = Float32[], Float32[]
    for y in 0:(NYEARS - 1), step in SAMPLE_STEPS
        f = forcing(y * N + step, s.start_year + y, r)
        push!(co2, f.CO2)
        push!(solar, f.sw_solar_forcing)
    end
    return (co2_ctrl = ini.CO2_ctrl, start_year = s.start_year, output = s.output,
            co2_table = co2_table, solar_table = _solar_name(s.solar),
            solar_row = s.solar isa SolarTable ? s.solar.index : nothing, boundary = _boundary_name(s.surface),
            processes = _parts(config.processes), hydrology = _parts(config.hydrology),
            corrections = repr(config.corrections),
            static_mask = rle(vec(static_mask)), dynamic_mask = rle(vec(fields.co2_part)),
            co2 = rle(co2), solar = rle(solar))
end
