# =============================================================================
# parameter_sweep.jl - CO2-concentration sensitivity sweep.
#
# Control at 280 ppm (the control CO2 of :custom_co2); scenario phase run at
# each level in co2_grid, the levels side by side with `run_ensemble`. Records
# scenario-minus-control anomalies of Ts, ice extent and precipitation as
# area-weighted global means.
#
# Run as script:  julia --project=. -t 4,0 examples/parameter_sweep.jl [data_dir]
# Or from REPL:   include("examples/parameter_sweep.jl"); parameter_sweep("data_dir")
#
# JLD2 data is not committed; pass its directory as argument or via GREB_DATA.
# Start Julia with several threads (`-t 4,0`) for the levels to run in parallel;
# each running level holds a copy of the input fields, about 370 MB.
# =============================================================================

using GREBClimate
using Statistics

"""
    parameter_sweep(jld2_dir; co2_grid=default_co2_grid(), spinup=3, ctrl=5, scnr=100)

Run the :custom_co2 experiment at each CO2 level (ppm) in `co2_grid`, with the
control at that preset's 280 ppm. Each level has its own CO2 table and runs on
its own copy of the fields. Returns a Vector of
(co2, Ts_anom, ice_anom, precip_anom) NamedTuples and writes
examples/parameter_sweep_results.csv.

`result.scnr` holds one `MonthlyRecord` per month and is already an anomaly
against the control's final-year monthly climatology, so no further subtraction
happens here. Each anomaly below is the mean over the scenario's final 12
records, i.e. its final year.
"""
function parameter_sweep(jld2_dir::AbstractString;
                          co2_grid::AbstractVector{<:Real}=default_co2_grid(),
                          spinup::Int=3, ctrl::Int=5, scnr::Int=100)

    if !isdir(jld2_dir)
        @warn """
        JLD2 data directory not found: $jld2_dir
        Pass it as an argument:  julia --project=. examples/parameter_sweep.jl <dir>
        or set the GREB_DATA environment variable. See DATA_README.md.
        """
        return nothing
    end

    scnr >= 12 || throw(ArgumentError("scnr must be >= 12 to take a final-year mean, got $scnr"))

    println("Loading GREB dataset from: ", jld2_dir)
    fields = load_climatology(jld2_dir; dataset=:ncep)
    run = RunSpec(ctrl=ctrl, scnr=scnr)

    # :custom_co2's scenario clock starts at 1950 and advances one year per
    # simulated year, so each table needs an entry per scenario year.
    years = 1950:(1950 + scnr - 1)

    # The final-year mean of an anomaly field, over the scenario's last 12 records
    final_year(result, field) = mean(global_mean(getproperty(rec, field)) for rec in @view result.scnr[end-11:end])

    anomalies = mktempdir() do tmpdir
        configs = map(enumerate(co2_grid)) do (i, co2)
            co2_path = joinpath(tmpdir, "co2_$(i).txt")
            open(co2_path, "w") do io
                for yr in years
                    println(io, yr, " ", co2)
                end
            end
            preset(:custom_co2; path=co2_path, corrections=SpinUp(spinup))
        end
        println("Running $(length(configs)) CO2 levels (spinup=$spinup, ctrl=$ctrl, scnr=$scnr years), ",
                "$(min(Threads.nthreads(), length(configs))) at a time...")
        run_ensemble(run, configs; fields, jld2_dir) do result
            (Ts_anom=final_year(result, :Ts), ice_anom=final_year(result, :ice),
             precip_anom=final_year(result, :precip))
        end
    end

    results = [(co2=co2, a...) for (co2, a) in zip(co2_grid, anomalies)]

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
