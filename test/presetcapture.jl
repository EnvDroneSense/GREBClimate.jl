# What each experiment preset imposes, captured without the dataset. Used by
# test_presets.jl and by tools/validation/preset_reference.jl, which writes the
# reference in test/data/preset_reference.jl. Needs testutils.jl.
#
# Start year, output mode, the tables a preset loads and its boundary forcing
# are decided inline in greb_model!; they are restated here from that code.

const PRESETS = [sort!(collect(keys(GREBClimate._EXPERIMENT_OVERRIDES))); :decon_mean_climate; :decon_2xco2]
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

function capture_preset(p::Symbol)
    cfg = create_experiment_config(p)
    fields = synthetic_fields()
    ini = quiet() do
        init_model!(cfg, fields)
    end
    static_mask = copy(fields.co2_part)
    ice = zeros(Float32, X, Y, 12)
    ice[:, abs.(GREBClimate.lat_grid) .> 60, :] .= 1.0f0
    GREBClimate.apply_dynamic_co2_mask!(cfg, fields, ice)

    orbital = p in (:obliquity, :eccentricity, :earth_sun_distance)
    start_year = orbital ? 1 : p === :historical_co2 ? 1850 : 1950
    co2_table = p === :custom_co2 ? :custom :
                p in GREBClimate._CO2_SCENARIO_SYMBOLS ? get(GREBClimate._CO2_SCENARIO_KEY, p, p) : :none
    if co2_table !== :none   # stand-in for the table greb_model! loads
        cfg.co2_scenario = Dict(y => 280.0f0 + 0.5f0 * (y - 1850) for y in start_year:(start_year + NYEARS))
    end

    co2, solar = Float32[], Float32[]
    for y in 0:(NYEARS - 1), s in SAMPLE_STEPS
        f = forcing(y * N + s, start_year + y, cfg, fields, ice)
        push!(co2, f.CO2)
        push!(solar, f.sw_solar_forcing)
    end
    return (co2_ctrl = ini.CO2_ctrl, start_year = start_year,
            output = orbital ? :absolute : :anomaly,
            co2_table = co2_table,
            solar_table = get(GREBClimate._SOLAR_SWAP_FORCING_TYPE, p, :none),
            boundary = p in (:rcp85, :elnino, :lanina, :sst_plus1) ? p : :none,
            static_mask = rle(vec(static_mask)), dynamic_mask = rle(vec(fields.co2_part)),
            co2 = rle(co2), solar = rle(solar))
end
