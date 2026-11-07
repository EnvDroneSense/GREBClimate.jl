# Which processes carry the warming of a doubled CO2?
#
#   julia --project=. -t 4,0 examples/feedback_deconstruction.jl [data_dir]
#
# Runs CO2 doubled once with the full model and once per process switched off,
# all side by side with `run_ensemble`, and prints the warming of the last
# scenario year for each. The table is also written to
# examples/out/feedback_deconstruction.csv. Members with different physics get
# their own spin-up. The differences do not add up: processes interact.
# Takes about a minute with 4 threads; each member holds 370 MB of fields.

using GREBClimate

const CTRL, SCNR = 3, 30    # years

# name => process switches that differ from the full model
const MEMBERS = [
    "full model" => (;),
    "no ice albedo" => (ice_albedo=false,),
    "uniform clouds" => (clouds=:uniform,),
    "no moisture convergence" => (moisture_convergence=false,),
    "no deep ocean" => (ocean=:mixed_layer,),
    "no transport" => (transport=false,),
    "no water cycle" => (hydrology=:none,)
]

function main(args)
    dir = greb_data_dir(isempty(args) ? nothing : args[1])
    fields = load_climatology(dir; dataset=:ncep)
    configs = [preset(:co2_double; processes=switches) for (_, switches) in MEMBERS]

    # The reducer runs inside each task: only the warming and a finite check are kept
    warming = run_ensemble(RunSpec(ctrl=CTRL, scnr=SCNR), configs; fields, jld2_dir=dir) do result
        last_year = result.scnr[(end-11):end]
        (K=sum(global_mean(r.Ts) for r in last_year) / 12,
            finite=all(r -> all(isfinite, r.Ts), result.scnr))
    end

    base = warming[1].K
    mkpath(joinpath(@__DIR__, "out"))
    csv = joinpath(@__DIR__, "out", "feedback_deconstruction.csv")
    open(csv, "w") do io
        println(io, "member,warming_K,difference_K")
        for ((name, _), w) in zip(MEMBERS, warming)
            println(io, name, ",", w.K, ",", w.K - base)
        end
    end

    println(rpad("member", 26), lpad("warming [K]", 12), lpad("vs full [K]", 13))
    for ((name, _), w) in zip(MEMBERS, warming)
        flag = w.finite ? "" : "  (run diverged)"
        println(rpad(name, 26), lpad(round(w.K, digits=2), 12), lpad(round(w.K - base, digits=2), 13), flag)
    end
    println("\nwrote ", csv)
    return warming
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
