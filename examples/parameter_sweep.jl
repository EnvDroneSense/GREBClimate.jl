# =============================================================================
# parameter_sweep.jl - CO2-concentration sensitivity sweep.
#
# Control held at the model's 340 ppm baseline; scenario phase run at each
# level in co2_grid via :custom_co2. Records scenario-minus-control anomalies
# in global-mean Ts, ice extent and precipitation.
#
# Run as script:  julia --project=. examples/parameter_sweep.jl [data_dir]
# Or from REPL:   include("examples/parameter_sweep.jl"); parameter_sweep("data_dir")
#
# JLD2 data is not committed; pass its directory as argument or via GREB_DATA.
# =============================================================================

using GREBClimate
using Statistics

"""
    parameter_sweep(jld2_dir; co2_grid=default_co2_grid(), flux=3, ctrl=5, scnr=100)

Run the :custom_co2 experiment at each CO2 level (ppm) in `co2_grid`, with the
control held at the model's 340 ppm baseline. Each grid point gets its own
config/fields instance and CO2 table. Returns a Vector of
(co2, Ts_anom, ice_anom, precip_anom) NamedTuples and writes
examples/parameter_sweep_results.csv.

`result.scnr` holds one `MonthlyRecord` per month and is already an anomaly
against the control's final-year monthly climatology (`apply_scenario_anomalies`
in src/postprocess.jl), so no further subtraction happens here. Each anomaly
below is the mean over the scenario's final 12 records, i.e. its final year.
"""
function parameter_sweep(jld2_dir::AbstractString;
                          co2_grid::AbstractVector{<:Real}=default_co2_grid(),
                          flux::Int=3, ctrl::Int=5, scnr::Int=100)

    if !isdir(jld2_dir)
        @warn """
        JLD2 data directory not found: $jld2_dir
        Pass it as an argument:  julia --project=. examples/parameter_sweep.jl <dir>
        or set the GREB_DATA environment variable. See DATA_README.md.
        """
        return nothing
    end

    println("Loading GREB dataset from: ", jld2_dir)
    fields_template = load_greb_jld2!(jld2_dir; dataset=:ncep)

    scnr >= 12 || throw(ArgumentError("scnr must be >= 12 to take a final-year mean, got $scnr"))

    run = RunSpec(flux=flux, ctrl=ctrl, scnr=scnr)
    results = NamedTuple{(:co2, :Ts_anom, :ice_anom, :precip_anom),
                          Tuple{Float64,Float64,Float64,Float64}}[]

    # :custom_co2's scenario clock starts at 1950 (src/model.jl:426) and
    # advances 1 year per step, so the table needs an entry per scenario year.
    years = 1950:(1950 + scnr - 1)

    mktempdir() do tmpdir
        for (i, co2) in enumerate(co2_grid)
            println("[$i/$(length(co2_grid))] CO2 = $(round(co2, digits=1)) ppm ",
                    "(flux=$flux, ctrl=$ctrl, scnr=$scnr years)...")

            co2_path = joinpath(tmpdir, "co2_$(i).txt")
            open(co2_path, "w") do io
                for yr in years
                    println(io, yr, " ", co2)
                end
            end

            cfg = create_experiment_config(:custom_co2; co2_path=co2_path)
            fields = deepcopy(fields_template)   # each grid point mutates its own state

            try
                result = greb_model!(run, cfg; jld2_dir=jld2_dir, fields=fields)

                scnr_final_year = @view result.scnr[end-11:end]  # last 12 monthly records
                Ts_anom = mean(mean(rec.Ts) for rec in scnr_final_year)
                ice_anom = mean(mean(rec.ice) for rec in scnr_final_year)
                precip_anom = mean(mean(rec.precip) for rec in scnr_final_year)

                push!(results, (co2=co2, Ts_anom=Ts_anom,
                                 ice_anom=ice_anom, precip_anom=precip_anom))
            catch err
                @error "sweep point failed" co2 exception=(err, catch_backtrace())
                push!(results, (co2=co2, Ts_anom=NaN, ice_anom=NaN, precip_anom=NaN))
            end
        end
    end

    out_path = joinpath(@__DIR__, "parameter_sweep_results.csv")
    open(out_path, "w") do io
        println(io, "co2_ppm,Ts_anom_K,ice_anom_frac,precip_anom")
        for r in results
            println(io, "$(r.co2),$(r.Ts_anom),$(r.ice_anom),$(r.precip_anom)")
        end
    end
    println("\nSaved ", length(results), " rows to ", out_path)

    return results
end

"""
    default_co2_grid()

20 log-spaced points spanning exactly the model's :co2_half (170 ppm) to
:co2_10x (3400 ppm) presets.
"""
default_co2_grid() = exp10.(range(log10(170.0), log10(3400.0), length=20))

if abspath(PROGRAM_FILE) == @__FILE__
    # Directory passed as first CLI arg, else GREB_DATA, else package default.
    jld2_dir = greb_data_dir(isempty(ARGS) ? nothing : ARGS[1])
    parameter_sweep(jld2_dir)
end
