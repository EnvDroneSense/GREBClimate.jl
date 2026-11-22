# Timing, comparison and profiling of GREBClimate.jl: the command line of Bench.jl.
#
#   julia --project=benchmark -t 2,0 benchmark/run.jl <mode> [--flag=value ...]
#
#   time      time control (+ scenario) runs: --preset=NAME (full_model), --ctrl=N (1),
#             --scnr=N (0), --reps=N (3), --out=FILE (run times, one per line)
#   threads   the same run at -t 1, 2, 3, 4 in subprocesses; --threads=1,2 sets the list
#   kernels   time each physics stage of one step: --seconds=S per stage (0.5)
#   compare   circulation!(Ta) against circulation!(q), interleaved, with a self-check:
#             --seconds=S per variant (5); --counters adds instructions and cycles per
#             call (Linux, needs LinuxPerf in your own environment). Exits 1 when the self-check fails
#   profile   sample control runs: --years=N (25), --samples=N (10000), --out=DIR,
#             --viewer=pprof; run at -t 1
#   allocs    record every allocation of one run: --years=N (1), --out=DIR
#   dispatch  check one model step for runtime dispatch with JET; static; --out=DIR
#   diff      the change in shares between two profile directories: diff <before> <after>
#
# Every mode but `diff` takes --data=DIR (default: the local dataset, never downloaded).
# The old names year, years, stages, step and setup still work for now.
# For a comparison of your own code, call Bench.compare_kernels from a script.

include(joinpath(@__DIR__, "..", "tools", "ToolBox.jl"))
using .ToolBox
include(joinpath(@__DIR__, "Bench.jl"))

const TIMING = ("preset", "ctrl", "scnr", "reps")
const FLAGS = Dict(
    "time" => (TIMING..., "out", "data"),
    "threads" => (TIMING..., "threads", "data"),
    "kernels" => ("seconds", "data"),
    "compare" => ("seconds", "counters", "data"),
    "profile" => ("years", "samples", "out", "viewer", "data"),
    "allocs" => ("years", "out", "data"),
    "dispatch" => ("out", "data"),
    "diff" => ()
)

# Old mode name => (mode, flags it implies). The flags come first, so the command line overrides them.
const ALIASES = Dict(
    "year" => ("time", String[]),
    "years" => ("time", ["--ctrl=10", "--scnr=10", "--reps=1"]),
    "stages" => ("kernels", String[]),
    "step" => ("profile", String[]),
    "setup" => ("profile", ["--years=1"])
)

function expand_alias(args::Vector{String})
    (isempty(args) || !haskey(ALIASES, args[1])) && return args
    mode, implied = ALIASES[args[1]]
    return [mode; implied; args[2:end]]
end

"Run the mode named by `args`; returns the exit code."
function main(args::Vector{String})
    mode, o, positional = parse_args(expand_alias(args), FLAGS, "time")
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
    out = get(o, "out", nothing)

    if mode == "time"
        Bench.time_years(dir; experiment, ctrl, scnr, reps, out)
    elseif mode == "threads"
        counts = Tuple(parse.(Int, split(get(o, "threads", "1,2,3,4"), ',')))
        Bench.sweep_threads(dir; thread_counts=counts, experiment, ctrl, scnr, reps)
    elseif mode == "kernels"
        Bench.time_stages(dir; experiment, seconds=floatflag(o, "seconds", 0.5))
    elseif mode == "compare"
        Bench.compare_circulation(dir; seconds=floatflag(o, "seconds", 5), counters=haskey(o, "counters")) || return 1
    elseif mode == "profile"
        Bench.profile_cpu(dir; experiment, years=intflag(o, "years", 25; min=1),
            target=intflag(o, "samples", Bench.TARGET_SAMPLES; min=1), out_dir=out, viewer=get(o, "viewer", nothing))
    elseif mode == "allocs"
        Bench.profile_allocs(dir; experiment, years=intflag(o, "years", 1; min=1), out_dir=out)
    else  # dispatch
        Bench.profile_dispatch(dir; out_dir=out)
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
