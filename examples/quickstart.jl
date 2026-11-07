# Quick start: run the model once and read the result.
#
#   julia --project=. -t 2,0 examples/quickstart.jl [data_dir]
#   julia --project=. examples/quickstart.jl --data-free
#
# Runs a doubling of CO2 for 10 years after a 2-year control, prints the warming
# per year, then runs it again with the ice-albedo feedback switched off. Without
# a data_dir it uses GREB_DATA, a local greb_input_data/ or a download (asks
# first). `--data-free` runs without the dataset to check that the calls work;
# its numbers mean nothing.

using GREBClimate

# Annual mean of the global mean of `field`, one value per year of monthly records
function yearly(records, field)
    return [sum(global_mean(getproperty(r, field)) for r in records[(12y-11):(12y)]) / 12
            for y in 1:(length(records)÷12)]
end

function main(args)
    if "--data-free" in args
        cfg = preset(:co2_double; corrections=NoCorrections())
        result = greb_model!(RunSpec(ctrl=1, scnr=1), cfg; jld2_dir="", allow_uninitialized=true)
        println("data-free run: ", length(result.ctrl), " control and ", length(result.scnr),
            " scenario months, fields of size ", size(result.ctrl[1].Ts))
        return result
    end

    # 1. The input data: a directory, loaded into `fields`
    dir = greb_data_dir(isempty(args) ? nothing : args[1])
    fields = load_climatology(dir; dataset=:ncep)

    # 2. The experiment: a named preset, here CO2 doubled over the control
    cfg = preset(:co2_double)

    # 3. The run: 2 control years, then 10 scenario years (after a 3-year spin-up)
    run = RunSpec(ctrl=2, scnr=10)
    result = greb_model!(run, cfg; jld2_dir=dir, fields=fields)

    # 4. The result: one record per month, each with 96x48 fields. The scenario
    #    records are anomalies against the control's last year.
    println("control: global-mean Ts ", round(sum(yearly(result.ctrl, :Ts)) / 2, digits=1), " K")

    # 5. A change to the experiment: the same run with ice albedo switched off
    no_ice = greb_model!(run, preset(:co2_double; processes=(ice_albedo=false,));
        jld2_dir=dir, fields=fields)

    warming, warming_no_ice = yearly(result.scnr, :Ts), yearly(no_ice.scnr, :Ts)
    println("\nyear   warming [K]   without ice albedo [K]")
    for y in eachindex(warming)
        println(lpad(y, 4), lpad(round(warming[y], digits=2), 14), lpad(round(warming_no_ice[y], digits=2), 18))
    end
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
