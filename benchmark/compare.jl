# Interleaved timing of code variants in one process. Uses Chairmarks.@b only.
# Needs the benchmark environment (julia --project=benchmark) and `common.jl`
# included first.

using Chairmarks
using Printf
using Random
using Statistics

const ROUNDS = 20             # timed chunks per variant
const SELF_CHECK_LIMIT = 0.01 # the baseline copy's median ratio must be within 1 % of 1

"Fastest time in seconds of repeated calls to `f` over about `seconds`."
chunk_time(f, seconds::Real) = (@b f seconds=seconds).time

"Why hardware counters cannot be read here, or `nothing` when they can."
function counters_unavailable()
    Sys.islinux() || return "counters need Linux"
    paranoid = parse(Int, strip(read("/proc/sys/kernel/perf_event_paranoid", String)))
    paranoid <= 2 || return "perf_event_paranoid is $paranoid; counters need 2 or lower"
    Base.find_package("LinuxPerf") === nothing && return "LinuxPerf is not installed; run Pkg.add(\"LinuxPerf\")"
    return nothing
end

"""
Compare `"name" => f` variants. The first pair is the baseline; a second copy
of it is added as `"<name> (copy)"`. Each of `rounds` rounds times every
variant for `seconds / rounds`, in a freshly shuffled order, so no variant
always runs last. A sample is the fastest call of one round's chunk.

Prints min and median per variant and the ratios to the baseline, appends the
records to `benchmark/results/`, and returns whether the copy's median ratio
is within 1 % of 1.000. When it is not, the run is not trustworthy.

With `counters=true` it also prints instructions and cycles per call (Linux
only, see `counters.jl`); where they cannot be read it says why and goes on.
"""
function compare_kernels(pairs::Pair{String}...; seconds::Real=5, rounds::Int=ROUNDS, counters::Bool=false)
    names = String[first(p) for p in pairs]
    fs = Any[last(p) for p in pairs]
    insert!(names, 2, names[1] * " (copy)")
    insert!(fs, 2, fs[1])
    for f in fs
        f()   # compiled before it is timed
    end

    calibration = machine_header()
    times = [Float64[] for _ in fs]
    for _ in 1:rounds, i in shuffle(eachindex(fs))
        push!(times[i], chunk_time(fs[i], seconds / rounds))
    end
    quiet = machine_footer(calibration)

    mins = minimum.(times)
    meds = median.(times)
    selfcheck = meds[2] / meds[1]
    ok = abs(selfcheck - 1) < SELF_CHECK_LIMIT
    @printf("%-26s %10s %8s %10s %8s\n", "variant", "min (us)", "x base", "med (us)", "x base")
    for i in eachindex(fs)
        @printf("%-26s %10.3f %8.3f %10.3f %8.3f\n", names[i], 1e6 * mins[i], mins[i] / mins[1], 1e6 * meds[i], meds[i] / meds[1])
        record_result("compare", names[i], times[i], quiet; ratios=(mins[i] / mins[1], meds[i] / meds[1]),
            selfcheck=i == 2 ? selfcheck : nothing)
    end
    println("self-check (copy / baseline, median; min gives ", round(mins[2] / mins[1], digits=3), "): ",
        round(selfcheck, digits=3), ok ? "  OK" : "  NOT OK - do not quote these ratios", "  [", rounds, " samples each]")
    counters && print_counters(names, fs)
    return ok
end

"Print instructions and cycles per call of every variant, refusing values whose counter did not run throughout."
function print_counters(names, fs)
    reason = counters_unavailable()
    if reason !== nothing
        println("counters: ", reason)
        return nothing
    end
    counter_values = include(joinpath(@__DIR__, "counters.jl"))
    values = [Base.invokelatest(counter_values, f) for f in fs]
    @printf("%-26s %12s %8s %12s %8s\n", "variant", "instr/call", "x base", "cycles/call", "x base")
    for i in eachindex(fs)
        instructions, cycles, fill = values[i]
        if fill < 1
            @printf("%-26s refused: fill rate %.3f is below 1.000 (pin to a P-core)\n", names[i], fill)
        else
            @printf("%-26s %12.0f %8.4f %12.0f %8.4f\n", names[i], instructions, instructions / values[1][1], cycles, cycles / values[1][2])
        end
    end
    return nothing
end
