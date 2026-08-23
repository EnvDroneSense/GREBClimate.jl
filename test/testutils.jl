# Shared fixtures and helpers for the test suite.
#
# `allow_download=false` so running the tests can never trigger the 353 MB
# DataDeps download; data-dependent testsets @test_skip when this is absent.
const DATA_DIR = something(greb_data_dir(; allow_download = false),
                           joinpath(@__DIR__, "..", "greb_input_data"))

const X, Y, N = GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr

"""Run `f` with stdout muted - the model prints a progress line per year."""
quiet(f) = redirect_stdout(devnull) do
    f()
end

"""Run `f(dir)` in a fresh temp directory, removed afterwards even on failure."""
function with_tempdir(f)
    dir = mktempdir()
    try
        f(dir)
    finally
        rm(dir; recursive = true, force = true)
    end
end

"""
Write a minimal `scenario/ipcc_scenarios.jld2` into `dir`. `table` maps the
on-disk scenario key ("rcp6", not "rcp60") to a year=>ppm Dict.
"""
function write_ipcc_scenarios(dir, table)
    mkpath(joinpath(dir, "scenario"))
    GREBClimate.jldopen(joinpath(dir, "scenario", "ipcc_scenarios.jld2"), "w") do file
        file["scenarios"] = table
    end
    return dir
end

"""
Write the three `solar_scenarios/*.jld2` tables the paleo/orbital experiments
swap in at scenario start, filled with a recognisable sentinel.
"""
function write_solar_scenarios(dir; sentinel = 999.0)
    mkpath(joinpath(dir, "solar_scenarios"))
    GREBClimate.jldopen(joinpath(dir, "solar_scenarios", "solar_paleo.jld2"), "w") do file
        file["data"] = fill(sentinel, Y, N)
        file["dim_names"] = ["lat", "time"]
    end
    for which in ("obliquity", "eccentricity")
        GREBClimate.jldopen(joinpath(dir, "solar_scenarios", "solar_" * which * ".jld2"), "w") do file
            file["data"] = fill(sentinel, 1, Y, N)
            file["dim_names"] = ["index", "lat", "time"]
            file["coords"] = Dict(1 => [0.0])
        end
    end
    return dir
end

"""
A `ClimateFields` non-degenerate enough to exercise the full physics pipeline
without the real dataset: finite mixed-layer depth, land in the western half,
and non-zero winds/humidity so circulation actually transports something.
"""
function synthetic_fields()
    f = ClimateFields()
    f.mldclim .= 50.0f0
    f.z_topo[1:(X - 48), :] .= 100.0f0
    for k in 1:N, j in 1:Y, i in 1:X
        f.Tclim[i, j, k] = 288.0f0 - 40.0f0 * abs(j - Y / 2) / (Y / 2)
        f.uclim[i, j, k] = 5.0f0 * sinpi(2 * j / Y)
        f.vclim[i, j, k] = 2.0f0 * cospi(2 * i / X)
        f.qclim[i, j, k] = 0.005f0
        f.wsclim[i, j, k] = 6.0f0
    end
    f.Toclim .= 283.0f0
    f.cldclim .= 0.5f0
    f.swetclim .= 0.4f0
    f.uclim_p .= max.(f.uclim, 0.0f0)
    f.uclim_m .= min.(f.uclim, 0.0f0)
    f.vclim_p .= max.(f.vclim, 0.0f0)
    f.vclim_m .= min.(f.vclim, 0.0f0)
    return f
end

gmean(x) = sum(x) / length(x)
