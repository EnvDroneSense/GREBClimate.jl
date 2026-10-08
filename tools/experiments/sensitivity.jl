# Effective climate sensitivity of an abrupt-CO2 run, from its monthly records.
#
#   julia --project=. -t 2,0 tools/experiments/sensitivity.jl [preset] [scenario years] [spin-up years] [control years]
#
# Runs `preset` (default co2_double) for the given scenario years (default 100),
# on a spin-up of 3 years and a control of 5, with the NCEP climatology, and
# fits the top-of-atmosphere net-flux change (absorbed shortwave minus `olr`)
# against the surface-temperature change over the annual global means (Gregory
# regression). Prints the forcing, the feedback, the effective sensitivity and
# r2, the warming the run reached and the net flux still left at the end.
#
# The regression value can sit below the reached warming: `sw - olr` is not the
# whole budget of the scenario anomaly. For `co2_double`, 100 years on the NCEP
# data, the numbers are 5.06 W/m2, -1.88 W/m2/K, 2.69 K, reached 3.00 K.
# Needs the local dataset and a preset whose scenario output is an anomaly.

include(joinpath(@__DIR__, "helpers.jl"))
using GREBClimate
using Printf
using .ExperimentTools: summarize, gregory

function main(args)
    name = Symbol(get(args, 1, "co2_double"))
    years = parse(Int, get(args, 2, "100"))
    spinup = parse(Int, get(args, 3, "3"))
    ctrl = parse(Int, get(args, 4, "5"))

    dir = greb_data_dir(; allow_download=false)
    dir === nothing && error("no local dataset: this tool never downloads it")
    fields = load_climatology(dir; dataset=:ncep)
    result = greb_model!(RunSpec(ctrl=ctrl, scnr=years), preset(name; corrections=SpinUp(spinup));
                         jld2_dir=dir, fields=fields)
    result.scnr_anomaly || error("the scenario of `$name` is not an anomaly; the estimate needs the change against the control")

    s = summarize(result)
    for phase in (:ctrl, :scnr)
        flag = s.out_of_range[phase]
        flag === nothing || @warn "the $phase run left the allowed range" flag
    end
    g = gregory(s.scnr.Ts, s.scnr.net)
    @printf("%s, %d scenario years, spin-up %d, control %d\n", name, years, spinup, ctrl)
    @printf("  forcing                %7.2f W/m2\n", g.forcing)
    @printf("  feedback               %7.2f W/m2/K\n", g.feedback)
    @printf("  effective sensitivity  %7.2f K\n", g.sensitivity)
    @printf("  r2                     %7.3f\n", g.r2)
    @printf("  warming reached        %7.2f K\n", g.reached)
    @printf("  net flux remaining     %7.2f W/m2\n", g.remaining)
    return g
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
