# Hardware counters through LinuxPerf (Linux only). Included by `compare.jl` at
# run time, after it has checked that counters can be read.

using LinuxPerf

"""
Instructions and cycles per call of `f`, and the lowest fill rate seen. Each
of `repeats` measurements covers `calls` calls; the counts are medians. A fill
rate below 1 means a counter was not running the whole time (for example on an
E-core), so its counts are scaled estimates.
"""
function counter_values(f, calls::Int=200, repeats::Int=10)
    run_calls() =
        for _ in 1:calls
            f()
        end
    run_calls()
    instructions = Float64[]
    cycles = Float64[]
    fill = 1.0
    for _ in 1:repeats
        thread = (@pstats "(cpu-cycles,instructions)" run_calls()).threads[1]
        fill = min(fill, LinuxPerf.fillrate(thread["instructions"]), LinuxPerf.fillrate(thread["cpu-cycles"]))
        push!(instructions, LinuxPerf.scaledcount(thread["instructions"]) / calls)
        push!(cycles, LinuxPerf.scaledcount(thread["cpu-cycles"]) / calls)
    end
    return median(instructions), median(cycles), fill
end

counter_values   # the value of this file, for `include`
