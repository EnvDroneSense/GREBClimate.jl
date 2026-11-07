# Climate response over a range of CO2 levels.
#
#   julia --project=. -t 4,0 tools/experiments/co2_sweep.jl [levels] [scenario years] [spin-up years] [control years]
#
# Holds CO2 at each of `levels` values (default 5), log-spaced between 170 ppm
# (:co2_half) and 3400 ppm (:co2_10x), for the scenario years (default 30) after
# a spin-up (default 3) and a control at 280 ppm (default 5). Prints, and writes
# to tools/experiments/out/co2_sweep.csv, the anomaly of the last scenario year:
# global-mean Ts, ice cover and precipitation. For a full curve use
# `... 20 100`, which takes hours. Needs the local dataset and several threads;
# each running level holds 370 MB of fields.

include(joinpath(@__DIR__, "helpers.jl"))
using GREBClimate
using .ExperimentTools: annual_means

function main(args)
    levels = parse(Int, get(args, 1, "5"))
    years = parse(Int, get(args, 2, "30"))
    spinup = parse(Int, get(args, 3, "3"))
    ctrl = parse(Int, get(args, 4, "5"))
    years >= 12 || error("need at least 12 scenario years to take a last-year mean")

    dir = greb_data_dir(; allow_download=false)
    dir === nothing && error("no local dataset found")
    fields = load_climatology(dir; dataset=:ncep)
    co2 = levels == 1 ? [170.0] : exp10.(range(log10(170.0), log10(3400.0), length=levels))
    configs = [Config(scenario=Scenario(co2=ConstantCO2(c), control_co2=280), corrections=SpinUp(spinup))
               for c in co2]

    # Only the last-year anomalies leave each task
    last_year(result, field) = annual_means(result.scnr, result.scnr_time, field)[2][end]
    anomalies = run_ensemble(RunSpec(ctrl=ctrl, scnr=years), configs; fields, jld2_dir=dir) do result
        (Ts=last_year(result, :Ts), ice=last_year(result, :ice), precip=last_year(result, :precip))
    end

    mkpath(joinpath(@__DIR__, "out"))
    csv = joinpath(@__DIR__, "out", "co2_sweep.csv")
    open(csv, "w") do io
        println(io, "co2_ppm,Ts_anom_K,ice_anom_frac,precip_anom")
        for (c, a) in zip(co2, anomalies)
            println(io, c, ",", a.Ts, ",", a.ice, ",", a.precip)
        end
    end
    println(lpad("CO2 [ppm]", 10), lpad("Ts [K]", 10), lpad("ice", 10), lpad("precip", 10))
    for (c, a) in zip(co2, anomalies)
        println(lpad(round(c, digits=0), 10), lpad(round(a.Ts, digits=2), 10),
            lpad(round(a.ice, digits=4), 10), lpad(round(a.precip, digits=3), 10))
    end
    println("\nwrote ", csv)
    return anomalies
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
