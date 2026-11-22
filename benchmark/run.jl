# Timing, comparison and profiling of GREBClimate.jl: the command line of Bench.jl.
#
#   julia --project=benchmark -t 2,0 benchmark/run.jl <mode> [--flag=value ...]
#
#   time      time control (+ scenario) runs: --ctrl=N (1), --scnr=N (0), --reps=N (3),
#             --out=FILE (run times, one per line)
#   threads   the same run at -t 1, 2, 3, 4 in subprocesses; --threads=1,2 sets the list
#   kernels   time the kernels of one step (the kernel table of Bench.jl, `step!` is a
#             whole step): --only=hydro!,seaice! (all), --seconds=S per kernel (0.5)
#   compare   time two or more kernels of the table interleaved, with a self-check,
#             the first being the baseline: --only=A,B[,C], --seconds=S per variant (5);
#             --counters adds instructions and cycles per call (Linux, needs LinuxPerf
#             in your own environment). Exits 1 when the self-check fails
#   profile   sample control runs: --years=N (25), --samples=N (10000), --out=DIR; run at -t 1
#   allocs    record every allocation of one run: --years=N (1), --out=DIR
#   dispatch  check one model step for runtime dispatch with JET; static; --out=DIR
#   diff      the change in shares between two profile directories: diff <before> <after>
#
# Every mode but `diff` takes --preset=NAME (full_model) and --data=DIR (default: the
# local dataset, never downloaded). To time, profile or compare your own code, call
# Bench.compare, Bench.profile, Bench.allocs or Bench.dispatch with a closure.

include(joinpath(@__DIR__, "..", "tools", "ToolBox.jl"))
using .ToolBox
include(joinpath(@__DIR__, "Bench.jl"))

const FLAGS = Dict(
    "time" => ("preset", "data", "ctrl", "scnr", "reps", "out"),
    "threads" => ("preset", "data", "ctrl", "scnr", "reps", "threads"),
    "kernels" => ("preset", "data", "only", "seconds"),
    "compare" => ("preset", "data", "only", "seconds", "counters"),
    "profile" => ("preset", "data", "years", "samples", "out"),
    "allocs" => ("preset", "data", "years", "out"),
    "dispatch" => ("preset", "data", "out"),
    "diff" => ()
)

"The comma-separated list of flag `name`; empty when absent."
list(o, name) = haskey(o, name) ? String.(split(o[name], ',')) : String[]

"Run the mode named by `args`; returns the exit code."
function main(args::Vector{String})
    mode, o, positional = parse_args(args, FLAGS, "time")
    if mode == "diff"
        length(positional) == 2 || error("diff takes two report directories: the run before and the run after")
        Bench.compare_runs(positional[1], positional[2])
        return 0
    end
    isempty(positional) || error("the $mode mode takes no positional argument, got: ", join(positional, " "))

    dir = haskey(o, "data") ? o["data"] : data_dir()
    isdir(dir) || error("dataset directory not found: $dir")
    experiment = Symbol(get(o, "preset", "full_model"))
    ctrl, scnr, reps = intflag(o, "ctrl", 1; min=1), intflag(o, "scnr", 0), intflag(o, "reps", 3; min=1)
    years = intflag(o, "years", mode == "allocs" ? 1 : 25; min=1)
    out = get(o, "out", nothing)
    header = ["preset" => String(experiment), "years" => string(years)]

    if mode == "time"
        Bench.time_years(dir; experiment, ctrl, scnr, reps, out)
    elseif mode == "threads"
        counts = Tuple(parse.(Int, haskey(o, "threads") ? list(o, "threads") : ["1", "2", "3", "4"]))
        Bench.sweep_threads(dir; thread_counts=counts, experiment, ctrl, scnr, reps)
    elseif mode == "kernels"
        Bench.time_kernels(dir; experiment, only=list(o, "only"), seconds=floatflag(o, "seconds", 0.5))
    elseif mode == "compare"
        only = list(o, "only")
        length(only) >= 2 || error("compare needs --only=A,B with two or more kernels of: ",
            join(first.(Bench.kernel_table(dir, experiment)), ", "))
        variants = Bench.pick(Bench.kernel_table(dir, experiment), only)
        Bench.compare(variants...; seconds=floatflag(o, "seconds", 5), counters=haskey(o, "counters")) || return 1
    elseif mode == "profile"
        Bench.profile(Bench.workload(dir, experiment; ctrl=years); warmup=false, out_dir=out, info=header,
            description=Bench.workload_text(experiment, years), target=intflag(o, "samples", Bench.TARGET_SAMPLES; min=1),
            max_runs=max(8, 400 ÷ years))
    elseif mode == "allocs"
        Bench.allocs(Bench.workload(dir, experiment; ctrl=years); warmup=false, out_dir=out, info=header,
            description=Bench.workload_text(experiment, years))
    else  # dispatch
        step = Dict(Bench.kernel_table(dir, experiment))["step!"]
        Bench.dispatch(step; out_dir=out, info=["preset" => String(experiment)],
            description="one time_loop! step of $experiment")
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
