# Timing, comparison and profiling of GREBClimate.jl in one module. The command
# line is benchmark/run.jl; from the REPL:
#
#   include("benchmark/Bench.jl"); using .Bench
#
# Run it with `julia --project=benchmark` (Chairmarks). Wall-clock numbers of
# different sessions are not comparable; see benchmark/README.md.

# `Bench` needs the shared helpers of tools/ToolBox.jl.
isdefined(@__MODULE__, :ToolBox) || include(joinpath(@__DIR__, "..", "tools", "ToolBox.jl"))

module Bench

using GREBClimate
using GREBClimate: SWradiation!, LWradiation!, hydro!, deep_ocean!, circulation!, seaice!
using Chairmarks: @b
using Dates
using Printf
using Profile
using Random: shuffle
using Statistics: median
using ..ToolBox: quiet, repo_state

const REPO = normpath(joinpath(@__DIR__, ".."))
const JULIA_BIN = joinpath(Sys.BINDIR, Base.julia_exename())
const SRC_DIR = normpath(joinpath(REPO, "src"))
const PROFILES_DIR = joinpath(REPO, "benchmark", "profiles")

# ----- Machine state -----
#
# A timing only means something next to how fast the machine was at that
# moment. A fixed loop is timed before and after the benchmark and compared
# with the best time this machine has ever given for it.

const CALIBRATION_FILE = joinpath(REPO, "benchmark", "results", "calibration.txt")
const NOISE_LIMIT = 1.05   # a ratio above this marks the run as noisy

@noinline function calibration_work()
    s = 0.0
    for i in 1:20_000_000
        s += sqrt(Float64(i))
    end
    return s
end

"Seconds for the fixed calibration loop: the minimum of `n` runs."
function calibrate(n::Int=5)
    calibration_work()   # compiled before it is timed
    return minimum(@elapsed(calibration_work()) for _ in 1:n)
end

"The best calibration time recorded on this machine, or `nothing`."
function best_calibration()
    isfile(CALIBRATION_FILE) || return nothing
    return tryparse(Float64, strip(read(CALIBRATION_FILE, String)))
end

"Record `t` when it beats the best calibration time; returns the best."
function record_calibration(t::Float64)
    best = best_calibration()
    if best === nothing || t < best
        mkpath(dirname(CALIBRATION_FILE))
        write(CALIBRATION_FILE, string(t))
        return t
    end
    return best
end

"Processor load, power source and background processes, as text; `unknown` where the system does not say."
function system_state()
    if Sys.iswindows()
        script = raw"$l=(Get-CimInstance Win32_Processor | Measure-Object LoadPercentage -Average).Average; " *
                 raw"$b=(Get-CimInstance Win32_Battery).BatteryStatus; " *
                 raw"$p=@(Get-Process OneDrive,SearchIndexer -ErrorAction SilentlyContinue | ForEach-Object Name | Sort-Object -Unique) -join ','; " *
                 raw"Write-Output (($l, $b, $p) -join ';')"
        try
            load, battery, procs = split(readchomp(Cmd(["powershell", "-NoProfile", "-Command", script])), ';')
            power = isempty(battery) ? "mains (no battery)" : battery == "1" ? "battery" : "mains"
            return "processor load $(load) %, $power, background: $(isempty(procs) ? "none of OneDrive, SearchIndexer" : procs)"
        catch
            return "unknown"
        end
    end
    governor = try
        strip(read("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor", String))
    catch
        "unknown"
    end
    return "load average $(round(Sys.loadavg()[1], digits=2)) on $(Sys.CPU_THREADS) threads, " *
           "power profile $(power_profile()), governor $governor, cpu $(ccall(:sched_getcpu, Cint, ()))"
end

"The Linux power profile (`powerprofilesctl get`); `unknown` where the system does not say."
function power_profile()
    Sys.islinux() || return "unknown"
    try
        return readchomp(`powerprofilesctl get`)
    catch
        return "unknown"
    end
end


"Print the machine state before a benchmark; returns the calibration time."
function machine_header()
    t = calibrate()
    best = record_calibration(t)
    println("machine: ", system_state())
    println("calibration loop: ", round(1e3 * t, digits=2), " ms, ", round(t / best, digits=2),
        "x the best seen here (", round(1e3 * best, digits=2), " ms)")
    return t
end

"""
Print the machine state after a benchmark and whether the run was noisy:
the calibration loop against the best seen and against `before`, the power profile, and for
`times` the ratio of the median to the minimum.
"""
function machine_footer(before::Float64, times::Vector{Float64}=Float64[])
    after = calibrate()
    best = record_calibration(after)
    worst = max(before, after) / best
    reasons = String[]
    worst > NOISE_LIMIT && push!(reasons, "calibration $(round(worst, digits=2))x the best")
    power_profile() == "power-saver" && push!(reasons, "power profile power-saver")
    if length(times) >= 3
        spread = median(times) / minimum(times)
        println("median / minimum of the runs: ", round(spread, digits=2))
        spread > NOISE_LIMIT && push!(reasons, "runs spread $(round(spread, digits=2))x")
    end
    println("calibration loop after: ", round(1e3 * after, digits=2), " ms, ", round(after / best, digits=2), "x the best")
    println(isempty(reasons) ? "machine state: quiet" : "machine state: NOISY (" * join(reasons, "; ") * ") - do not quote these timings")
    return isempty(reasons)
end

# ----- Result records -----

const RESULTS_DIR = joinpath(REPO, "benchmark", "results")
const RECORD_COLUMNS = "date,commit,mode,name,threads,samples,min_s,median_s,ratio_min,ratio_median,selfcheck,state"

"""
Append one line to `benchmark/results/<date>.csv`, creating the file with its
header line when it is new. `ratios` is `(min, median)` against the baseline,
`nothing` for a single-variant mode; `selfcheck` is the baseline copy's median
ratio or `nothing`; `calm` is what `machine_footer` returned.
"""
function record_result(mode::AbstractString, name::AbstractString, times::Vector{Float64}, calm::Bool;
        ratios=nothing, selfcheck=nothing, threads::Int=Threads.nthreads())
    commit, dirty = repo_state()
    fields = (Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"), commit * (dirty ? "-dirty" : ""), mode, name,
        threads, length(times), minimum(times), median(times),
        ratios === nothing ? "" : ratios[1], ratios === nothing ? "" : ratios[2],
        selfcheck === nothing ? "" : selfcheck, calm ? "quiet" : "NOISY")
    mkpath(RESULTS_DIR)
    file = joinpath(RESULTS_DIR, Dates.format(today(), "yyyy-mm-dd") * ".csv")
    isfile(file) || write(file, RECORD_COLUMNS * "\n")
    open(io -> println(io, join(fields, ',')), file, "a")
    return nothing
end

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
function compare(pairs::Pair{String}...; seconds::Real=5, rounds::Int=ROUNDS, counters::Bool=false)
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
    calm = machine_footer(calibration)

    mins = minimum.(times)
    meds = median.(times)
    selfcheck = meds[2] / meds[1]
    ok = abs(selfcheck - 1) < SELF_CHECK_LIMIT
    @printf("%-26s %10s %8s %10s %8s\n", "variant", "min (us)", "x base", "med (us)", "x base")
    for i in eachindex(fs)
        @printf("%-26s %10.3f %8.3f %10.3f %8.3f\n", names[i], 1e6*mins[i], mins[i]/mins[1], 1e6*meds[i], meds[i]/meds[1])
        record_result("compare", names[i], times[i], calm; ratios=(mins[i] / mins[1], meds[i] / meds[1]),
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
            @printf("%-26s %12.0f %8.4f %12.0f %8.4f\n", names[i], instructions, instructions/values[1][1], cycles, cycles/values[1][2])
        end
    end
    return nothing
end

# ----- Workloads -----

"The climatology, workspace and slice-1 state that the kernel table calls the kernels with."
function stage_inputs(dir::AbstractString, cfg)
    fields = quiet(() -> load_climatology(dir; dataset=:ncep))
    r = resolve(cfg; jld2_dir=dir)
    CO2 = GREBClimate.init_model!(fields, r).CO2_ctrl
    slice = 1
    Ts = copy(fields.Ts_clim[:, :, slice])
    return (; fields, r, CO2, state=ModelState(), ws=ModelWorkspace(), slice, Ts, Ta=copy(Ts),
        To=copy(fields.To_clim[:, :, slice]), q=copy(fields.q_clim[:, :, slice]))
end

"Run `greb_model!` without its log output."
run_model(spec, cfg, dir, fields) = quiet(() -> greb_model!(spec, cfg; jld2_dir=dir, fields=fields))

"What a report says is run."
workload_text(experiment::Symbol, ctrl::Int=1, scnr::Int=0) = "$experiment, $ctrl control + $scnr scenario year(s), stored corrections"

"""
A closure that runs `ctrl` control years and `scnr` scenario years of the preset
`experiment` on the stored flux corrections (no spin-up). One warm-up year runs
first, so compilation is not timed. `greb_model!` restores `fields`, so every call
reuses them.
"""
function workload(dir::AbstractString, experiment::Symbol; ctrl::Int=1, scnr::Int=0)
    cfg = preset(experiment; corrections=Stored())
    fields = quiet(() -> load_climatology(dir; dataset=:ncep))
    spec = RunSpec(ctrl=ctrl, scnr=scnr)
    run_model(RunSpec(ctrl=1, scnr=0), cfg, dir, fields)
    GC.gc()
    return () -> run_model(spec, cfg, dir, fields)
end

"""
The kernels of one timestep as `name => call` on the slice-1 state of the preset
`experiment`, in the order the `kernels` mode prints them. `step!` is one whole
`time_loop!` step on its own copy of the state. `circulation!(q)` holds
`convergence!`. A name is what `--only` takes.
"""
function kernel_table(dir::AbstractString, experiment::Symbol=:full_model)
    (; fields, r, CO2, state, ws, slice, Ts, Ta, To, q) = stage_inputs(dir, preset(experiment))
    p = r.config.processes
    # `step!` advances its own copy, so it does not change what the kernels are called with.
    s = (; Ts=copy(Ts), Ta=copy(Ta), To=copy(To), q=copy(q), state=ModelState(), ws=ModelWorkspace(),
        acc=GREBClimate.MonthlyAccumulator(), records=GREBClimate.MonthlyRecord[],
        times=GREBClimate.RecordTime[], time=ModelTime(GREBClimate.control, 1970, 1))
    step = function ()
        empty!(s.records)   # a year of records would otherwise pile up while it is timed
        empty!(s.times)
        GREBClimate.time_loop!(s.time, CO2, s.Ts, s.Ta, s.To, s.q, s.records, s.times, fields, s.state, s.ws, s.acc, r)
    end
    return Pair{String,Function}[
        "circulation!(Ta)"=>()->circulation!(ws.dTa_crcl, Ta, GREBClimate.z_air, fields, ws, slice, p),
        "circulation!(q)"=>()->circulation!(ws.dq_crcl, q, GREBClimate.z_vapor, fields, ws, slice, p),
        "SWradiation!"=>()->SWradiation!(ws, Ts, fields, state, slice, p),
        "LWradiation!"=>()->LWradiation!(ws, Ts, Ta, q, CO2, fields, slice, p),
        "hydro!"=>()->hydro!(ws, Ts, q, fields, slice, p, r.hydrology),
        "deep_ocean!"=>()->deep_ocean!(ws, Ts, To, fields, slice, p),
        "seaice!"=>()->seaice!(fields, Ts, slice, p),
        "step!"=>step
    ]
end

"The entries of `table` named in `names`, in that order; all of them when `names` is empty."
function pick(table, names)
    isempty(names) && return table
    byname = Dict(table)
    for n in names
        haskey(byname, n) || error("unknown kernel $(repr(n)); expected one of ", join(first.(table), ", "))
    end
    return Pair{String,Function}[n => byname[n] for n in names]
end

# ----- Timing -----

"""
Time `reps` runs of `ctrl` control years and `scnr` scenario years of the preset
`experiment`; returns seconds per run. `out` is a file that receives the run
times, one per line.

`:full_model` holds CO2 constant across the whole scenario regardless of `scnr`;
for a forced change, pass a preset with a real CO2 path.
"""
function time_years(dir::AbstractString; experiment::Symbol=:full_model, ctrl::Int=1, scnr::Int=0, reps::Int=3,
        out=nothing)
    reps >= 1 || throw(ArgumentError("reps must be at least 1, got $reps"))
    ctrl >= 1 || throw(ArgumentError("ctrl must be at least 1, got $ctrl"))
    scnr >= 0 || throw(ArgumentError("scnr must be >= 0, got $scnr"))
    total_years = ctrl + scnr

    println("Threads.nthreads() = ", Threads.nthreads())
    println("ctrl=$ctrl scnr=$scnr ($total_years simulated years/rep), preset=$experiment")
    calibration = machine_header()
    run! = workload(dir, experiment; ctrl, scnr)

    times = Float64[]
    for r in 1:reps
        push!(times, @elapsed run!())
        println("  run $r: ", round(times[end], digits=3), " s  (",
            round(times[end] / total_years, digits=3), " s/simulated year)")
    end

    println("min: ", round(minimum(times), digits=3), " s  median: ", round(median(times), digits=3), " s  ",
        "= ", round(median(times) / total_years, digits=3), " s/simulated year (median) over ", total_years, " years")
    record_result("time", String(experiment), times, machine_footer(calibration, times))
    out === nothing || write(out, join(times, "\n"))
    return times
end

"`x` with `digits` decimals (trailing zeros kept), right-aligned in `width` characters."
function _fixed(x::Real, digits::Int, width::Int)
    whole, frac = split(string(round(x; digits)), '.')
    return lpad(whole * "." * rpad(frac, digits, '0'), width)
end

"""
Time the kernels of the table named in `only` (all when empty) for about `seconds`
each; returns `(name, fastest seconds per call)`. A kernel is timed in `ROUNDS`
chunks; the fastest call of each chunk is kept, and the minimum and median of those
are printed. With `step!` in the list, each row also shows its share of one step.
"""
function time_kernels(dir::AbstractString; experiment::Symbol=:full_model, only=String[], seconds::Real=0.5)
    seconds > 0 || throw(ArgumentError("seconds must be positive, got $seconds"))
    kernels = pick(kernel_table(dir, experiment), only)
    for (_, f) in kernels
        f()
    end
    calibration = machine_header()

    println("Kernel timing (about ", seconds, " s each, single workspace, ",
        Threads.nthreads(), " thread(s) available but unused here):")
    samples = [[chunk_time(f, seconds / ROUNDS) for _ in 1:ROUNDS] for (_, f) in kernels]
    results = [(name, minimum(t)) for ((name, _), t) in zip(kernels, samples)]

    step = findfirst(r -> first(r) == "step!", results)
    for ((name, per_call), t) in zip(results, samples)
        share = step === nothing ? "" : "   (" * _fixed(100 * per_call / last(results[step]), 1, 5) * "% of step!)"
        println("  ", rpad(name, 18), "min ", _fixed(per_call * 1e6, 2, 7), " µs   median ",
            _fixed(median(t) * 1e6, 2, 7), " µs", share)
    end
    println("  (convergence! is inside circulation!(q); step! also holds the update loop, output! and")
    println("   diagnostics!. For shares of a whole run use the profile mode)")
    calm = machine_footer(calibration)
    for ((name, _), t) in zip(kernels, samples)
        record_result("kernels", name, t, calm)
    end
    return results
end

"Time the `time` mode in a subprocess per thread count; returns seconds per run for each."
function sweep_threads(dir::AbstractString; thread_counts=(1, 2, 3, 4), experiment::Symbol=:full_model,
        ctrl::Int=1, scnr::Int=0, reps::Int=3)
    script = joinpath(@__DIR__, "run.jl")
    project = dirname(Base.active_project())

    results = Dict{Int,Vector{Float64}}()
    for n in thread_counts
        println("--- -t $n ---")
        out = tempname()
        run(`$JULIA_BIN --project=$project -t $n,0 $script time --data=$dir --preset=$experiment --ctrl=$ctrl --scnr=$scnr --reps=$reps --out=$out`)
        results[n] = parse.(Float64, readlines(out))
        rm(out)
    end

    base_n = first(thread_counts)
    base = median(results[base_n])
    println("\nSpeedup of the median vs -t $base_n:")
    for n in thread_counts
        r = results[n]
        println("  -t $n: min ", round(minimum(r), digits=3), " s  median ", round(median(r), digits=3),
            " s  (", round(base / median(r), digits=2), "x)")
    end
    return results
end

# ----- Profiling -----

const TARGET_SAMPLES = 10_000
const MIN_SAMPLES = 5_000        # below this only the largest rows mean anything
const PROFILE_BUFFER = 10_000_000
const META_WORDS = 6             # threadid, taskid, cycle clock, sleep state, 0, 0 close each sample
const CONSOLE_ROWS = 8           # rows of each table that a run prints

const NO_OWNER = "(outside the package)"

# ----- Frames and stacks -----

"One stack frame, reduced to what the reports need."
struct Frame
    func::String
    file::String     # with forward slashes
    line::Int
    from_c::Bool
    owned::Bool      # lies in the package's `src/`
end

function Frame(fr::Base.StackTraces.StackFrame)
    file = string(fr.file)
    Frame(clean_name(string(fr.func)), replace(file, '\\' => '/'), fr.line, fr.from_c,
        startswith(normpath(file), SRC_DIR))
end

"`name` without the `#name#12` wrapping of keyword bodies and closures, whose number changes between builds."
function clean_name(name::AbstractString)
    m = match(r"^#+([^#]+)#+\d*$", name)
    (m === nothing || all(isdigit, m.captures[1])) ? String(name) : String(m.captures[1])
end

"`fr` as it is named in the reports."
label(fr::Frame) = string(fr.func, " ", basename(fr.file), ":", fr.line)

"Add `by` to the count of `key`."
bump!(counts::AbstractDict, key, by::Int=1) = (counts[key]=get(counts, key, 0) + by; nothing)

"The entries of `counts`, largest first; equal counts in the order of their keys, so a report is reproducible."
ranked(counts::AbstractDict) = sort!(collect(counts); by=kv -> (-last(kv), first(kv)))

"The `top` largest entries of `counts` as `name (count)`, joined by commas."
function largest(counts::AbstractDict, top::Int=4)
    join((string(k, " (", v, ")") for (k, v) in first(ranked(counts), top)), ", ")
end

"Split profile data recorded with `include_meta=true` into `(frames, meta)` index ranges, one pair per sample."
function split_samples(data::Vector{UInt64})
    samples = Tuple{UnitRange{Int},UnitRange{Int}}[]
    start = 1
    for i in eachindex(data)
        # The sleep state before the two closing zeros is never zero (Profile's own rule).
        if data[i] == 0 && i - start + 1 >= META_WORDS && data[i-1] == 0 && data[i-2] != 0
            push!(samples, (start:(i-META_WORDS), (i-META_WORDS+1):i))
            start = i + 1
        end
    end
    return samples
end

"`data` with the frames outside the outermost package frame cut from every sample; `frames_at` maps an instruction pointer to its frames."
function trim_to_package(data::Vector{UInt64}, samples, frames_at::AbstractDict)
    owned = Set(ip for (ip, frames) in frames_at if any(fr -> fr.owned, frames))
    trimmed = UInt64[]
    sizehint!(trimmed, length(data))
    for (frames, meta) in samples
        # Frames are stored innermost first, so the outermost package frame is the last match.
        last_owned = findlast(i -> data[i] in owned, frames)
        keep = last_owned === nothing ? frames : first(frames):frames[last_owned]
        append!(trimmed, view(data, keep))
        append!(trimmed, view(data, meta))
    end
    return trimmed
end

"The distinct stacks as `frames => samples` pairs, innermost first, inlined frames included, cut at the outermost package frame."
function distinct_stacks(data::Vector{UInt64}, samples, frames_at::AbstractDict)
    counts = Dict{Vector{UInt64},Int}()
    for (frames, _) in samples
        bump!(counts, data[frames])
    end
    return map(collect(counts)) do (ips, count)
        stack = Frame[]
        for ip in ips
            append!(stack, frames_at[ip])
        end
        last_owned = findlast(fr -> fr.owned, stack)
        last_owned === nothing || resize!(stack, last_owned)
        stack => count
    end
end

"The innermost package frame of `stack` and its function, or `nothing` if there is none."
function owner(stack::Vector{Frame})
    i = findfirst(fr -> fr.owned, stack)
    i === nothing && return nothing
    # A loop under `@simd` or `@inbounds` shows as "macro expansion"; name it after the function around it.
    j = findnext(fr -> fr.owned && fr.func != "macro expansion", stack, i)
    return stack[i], (j === nothing ? stack[i].func : stack[j].func)
end

"The kind of cost the innermost frame `fr` of a sample stands for."
function classify(fr::Frame)
    f, file = fr.func, fr.file
    fr.owned && return "package code"
    occursin(r"VectorizationBase|LoopVectorization|SLEEFPirates", file) && return "vectorized kernel"
    (occursin(r"typeinf|type_infer|jl_compile|codegen|generate_fptr", f) || occursin(r"/[Cc]ompiler/", file)) &&
        return "compilation"
    occursin(r"apply_generic|jl_invoke|apply_iterate", f) && return "runtime dispatch"
    occursin(r"gc_|[Aa]lloc|Heap|^GenericMemory$|^Array$|^Memory$|^free$", f) && return "allocation and GC"
    occursin(r"^mem(cpy|move|set)$|copyto|^copy$|deepcopy|^fill!$", f) && return "copy and fill"
    (occursin(r"poptask|^wait$|task_done_hook|task_get_next|yieldto|^schedule$|jl_switch|NtWait|NtDelay", f) ||
     endswith(file, "/task.jl")) && return "task and wait"
    (occursin(r"jl_fs_|^uv_|^ios_|Nt(Read|Write|Create)File", f) || occursin(r"JLD2|iostream|filesystem|Mmap", file)) &&
        return "file I/O"
    (occursin(r"/special/|/math\.jl$|libm", file) || (fr.from_c && occursin(r"^(exp|log|pow|sin|cos|tan)", f))) &&
        return "scalar math"
    fr.from_c && return "runtime and system (C)"
    return "Base (array access, scalar ops)"
end

# ----- Markdown tables -----

"Write a padded Markdown table; the columns in `left` are left-aligned, the others right-aligned."
function md_table(io::IO, header::Vector{String}, rows; left=(1,))
    cell(s) = replace(string(s), "|" => "\\|")
    all_rows = [map(cell, header), (map(cell, r) for r in rows)...]
    widths = [maximum(length(r[j]) for r in all_rows) for j in eachindex(header)]
    line(r) = string("| ", join((j in left ? rpad(r[j], widths[j]) : lpad(r[j], widths[j]) for j in eachindex(r)), " | "), " |")
    println(io, line(all_rows[1]))
    println(io, "|", join((j in left ? string(":", "-"^(widths[j] + 1)) : string("-"^(widths[j] + 1), ":")
            for j in eachindex(header)), "|"), "|")
    for r in Iterators.drop(all_rows, 1)
        println(io, line(r))
    end
    println(io)
end

"`name` in backticks, so a Markdown viewer shows it as code."
code(name) = string("`", name, "`")

"The share of `count` in `n` samples and its standard error in percentage points, as strings."
function share(count::Int, n::Int)
    p = count / n
    return @sprintf("%.2f", 100p), @sprintf("%.2f", 100*sqrt(p * (1 - p) / n))
end

"`bytes` with a unit, for a table."
format_bytes(bytes::Integer) = bytes >= 2^20 ? @sprintf("%.2f MB", bytes/2^20) : @sprintf("%.2f kB", bytes/2^10)

# ----- Reports of a sampled run -----

"Count the samples of `stacks` per class, per package function (`self`: innermost package frame, `onstack`: anywhere) and per package line."
function tally(stacks)
    classes = Dict{String,Int}()
    class_frames = Dict{String,Dict{String,Int}}()
    class_owners = Dict{String,Dict{String,Int}}()
    self, onstack, lines = Dict{String,Int}(), Dict{String,Int}(), Dict{String,Int}()
    n = 0
    for (stack, count) in stacks
        isempty(stack) && continue
        n += count
        o = owner(stack)
        class = classify(first(stack))
        bump!(classes, class, count)
        bump!(get!(Dict{String,Int}, class_frames, class), label(first(stack)), count)
        bump!(get!(Dict{String,Int}, class_owners, class), o === nothing ? NO_OWNER : o[2], count)
        if o === nothing
            bump!(self, NO_OWNER, count)
            continue
        end
        fr, func = o
        bump!(self, func, count)
        bump!(lines, string(basename(fr.file), ":", fr.line, "  ", func), count)
        for name in unique!([f.func for f in stack if f.owned && f.func != "macro expansion"])
            bump!(onstack, name, count)
        end
    end
    return (; samples=n, classes, class_frames, class_owners, self, onstack, lines)
end

"Write `counts.tsv` (`table`, `row`, `value`), which `diff` reads."
function write_counts(path::AbstractString, t, commit::AbstractString)
    open(path, "w") do io
        println(io, "table\trow\tvalue")
        println(io, "run\tcommit\t", commit)
        println(io, "run\tsamples\t", t.samples)
        for (name, table) in (("class", t.classes), ("function", t.self), ("line", t.lines)), (row, count) in ranked(table)

            println(io, name, "\t", row, "\t", count)
        end
    end
end

"Write the class, function and line tables; `top` limits the rows (console), without it all rows of `mincount` samples or more."
function write_sampled(io::IO, t; mincount::Int, ms_per_run::Float64, top=nothing)
    n = t.samples
    rows_of(items) = top === nothing ? items : first(items, top)

    println(io, "## By class\n")
    println(io, "Each sample is classed by its innermost frame. Shares are of ", n, " samples; +/- is one standard error in percentage points.\n")
    md_table(io, ["Class", "Samples", "Share %", "+/-"],
        ([class, string(c), share(c, n)...] for (class, c) in ranked(t.classes)))
    if top === nothing
        println(io, "Largest frames of each class, and the package function they are charged to:\n")
        md_table(io, ["Class", "Largest frames", "Charged to"],
            ([class, largest(t.class_frames[class]), largest(t.class_owners[class])] for (class, _) in ranked(t.classes)))
    end

    println(io, "## By function\n")
    println(io, "Self: the function holds the innermost package frame. On stack: it is anywhere in the stack. ",
        "ms: the on-stack share of one run (", round(ms_per_run, digits=1), " ms).",
        top === nothing ? string(" Rows below ", mincount, " samples are left out but counted in the total.") : "", "\n")
    names = sort!(collect(union(keys(t.self), keys(t.onstack))); by=k -> (-get(t.self, k, 0), -get(t.onstack, k, 0), k))
    filter!(k -> max(get(t.self, k, 0), get(t.onstack, k, 0)) >= mincount, names)
    md_table(io, ["Function", "Self", "Share %", "+/-", "On stack", "Share %", "ms"],
        ([code(k), string(get(t.self, k, 0)), share(get(t.self, k, 0), n)...,
             string(get(t.onstack, k, 0)), share(get(t.onstack, k, 0), n)[1],
             @sprintf("%.2f", ms_per_run*get(t.onstack, k, 0)/n)] for k in rows_of(names));
        left=(1,))

    println(io, "## By line\n")
    println(io, "The innermost package frame. A `@turbo` loop is one row, at its `@turbo for` line.\n")
    lines = filter(kv -> last(kv) >= mincount, ranked(t.lines))
    md_table(io, ["Line", "Self", "Share %", "+/-"],
        ([code(line), string(c), share(c, n)...] for (line, c) in rows_of(lines)))
end

"Write `stacks.folded`, the collapsed-stack format flame graph tools read."
function write_folded(path::AbstractString, stacks)
    folded = Dict{String,Int}()
    for (stack, count) in stacks
        bump!(folded, join((label(fr) for fr in Iterators.reverse(stack)), ";"), count)
    end
    open(path, "w") do io
        for (key, count) in ranked(folded)
            println(io, key, " ", count)
        end
    end
end

"Write `Profile.print` output for `data` to `path` without truncating long lines."
function write_report(path::AbstractString, data, lidict; kwargs...)
    open(path, "w") do io
        Profile.print(IOContext(io, :displaysize => (200, 250)), data, lidict; kwargs...)
    end
end

# ----- Runs -----

"An optional package from the default environment (JET); not a dependency of this script."
function load_optional(name::AbstractString, uuid::AbstractString)
    try
        Base.require(Base.PkgId(Base.UUID(uuid), name))
    catch
        error("this needs $name.jl in the default environment: julia -e 'using Pkg; Pkg.add(\"$name\")'")
    end
end

"Make the report directory (`out_dir`, or a new numbered one under `benchmark/profiles/`); returns it and the `item => value` rows every report starts with."
function open_report(mode::AbstractString, description::AbstractString, out_dir)
    commit, dirty = repo_state()
    dir = out_dir
    if dir === nothing
        base = joinpath(PROFILES_DIR, string(Dates.today(), "-", commit, dirty ? "-dirty" : "", "-", mode))
        dir, n = base, 1
        while ispath(dir)
            n += 1
            dir = string(base, "-", n)
        end
    end
    mkpath(dir)
    info = ["mode" => "$mode ($description)",
        "commit" => commit * (dirty ? " (tracked files modified)" : ""),
        "date" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM"),
        "julia" => string(VERSION)]
    return dir, info
end

"Write `report.md` and print its top rows; `body(io, top)` writes the tables (`top === nothing`: all rows)."
function close_report(body, dir::AbstractString, title::AbstractString, info)
    head(io) = (println(io, "# ", title, "\n"); md_table(io, ["Item", "Value"], ([k, v] for (k, v) in info); left=(1, 2)))
    open(joinpath(dir, "report.md"), "w") do io
        head(io)
        body(io, nothing)
    end
    head(stdout)
    body(stdout, CONSOLE_ROWS)
    println("full report: ", joinpath(dir, "report.md"))
    return dir
end


"""
Sample calls of `f` with the profiler until `target` samples or `max_runs` calls;
one call is one run of the workload. Returns the report directory. `description`
and `info` (`item => value` pairs) go in the report header. `warmup` calls `f` once
unprofiled first, so compilation is not sampled; leave it off when `f` is warm.

Run at `-t 1`: on Windows the sampler records the first thread only.
"""
function profile(f; description::AbstractString="a closure", info=Pair{String,String}[],
        target::Int=TARGET_SAMPLES, max_runs::Int=50, out_dir=nothing, warmup::Bool=true)
    if Threads.nthreads() > 1
        @warn "Profiling with $(Threads.nthreads()) threads: on Windows only the first thread is sampled. Use -t 1."
    end
    warmup && f()
    GC.gc()

    Profile.clear()
    Profile.init(n=PROFILE_BUFFER, delay=0.001)
    elapsed = 0.0
    runs = 0
    while runs < max_runs &&
          (runs == 0 || length(split_samples(Profile.fetch(include_meta=true))) < target)
        elapsed += @elapsed @profile f()
        runs += 1
        Profile.is_buffer_full() && error("the profile buffer filled up; lower the run length or --samples")
    end

    data, lidict = Profile.retrieve(include_meta=true)
    samples = split_samples(data)
    n = length(samples)
    n > 0 || error("the profiler recorded no samples")
    threads_seen = sort!(unique(Int(data[first(meta)]) for (_, meta) in samples))
    ms_per_run = 1000 * elapsed / runs

    report, header = open_report("profile", description, out_dir)
    push!(header, info..., "threads" => "$(Threads.nthreads()) (sampled: $(join(threads_seen, ", ")))",
        "runs" => string(runs),
        "elapsed" => "$(round(elapsed, digits=2)) s ($(round(ms_per_run, digits=1)) ms per run)",
        "samples" => "$n ($(round(1000 * elapsed / n, digits=2)) ms per sample)")

    frames_at = Dict(ip => Frame.(frames) for (ip, frames) in lidict)
    trimmed = trim_to_package(data, samples, frames_at)
    stacks = distinct_stacks(data, samples, frames_at)
    t = tally(stacks)
    mincount = max(2, n ÷ 1000)   # rows below 0.1 percent are noise
    write_counts(joinpath(report, "counts.tsv"), t, last(header[2]))
    write_report(joinpath(report, "tree.txt"), trimmed, lidict; format=:tree, mincount, noisefloor=2)
    write_folded(joinpath(report, "stacks.folded"), stacks)

    close_report(report, "Profile", header) do io, top
        write_sampled(io, t; mincount, ms_per_run, top)
    end
    n >= MIN_SAMPLES ||
        @warn "Only $n samples (fewer than $MIN_SAMPLES): shares of small rows are not reliable. Raise --samples or --years."
    return report
end

"The allocations charged to one function or line."
mutable struct AllocRow
    count::Int
    bytes::Int
    types::Dict{String,Int}   # type name => bytes
end
AllocRow() = AllocRow(0, 0, Dict{String,Int}())

"Charge one allocation of `bytes` bytes of `type` to `row`."
function add!(row::AllocRow, type::AbstractString, bytes::Int)
    row.count += 1
    row.bytes += bytes
    bump!(row.types, type, bytes)
end


"""
Record every allocation of one call of `f`; returns the report directory.
`description`, `info` and `warmup` are as for `profile`.
"""
function allocs(f; description::AbstractString="a closure", info=Pair{String,String}[], out_dir=nothing,
        warmup::Bool=true)
    warmup && f()
    Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate=1 f()
    recorded = Profile.Allocs.fetch().allocs
    isempty(recorded) && error("the allocation profiler recorded nothing")

    frame_of = Dict{Base.StackTraces.StackFrame,Frame}()
    by_func = Dict{String,AllocRow}()
    by_line = Dict{String,AllocRow}()
    total, in_loop = AllocRow(), AllocRow()
    for a in recorded
        stack = [get!(() -> Frame(fr), frame_of, fr) for fr in a.stacktrace]
        type = first(string(a.type), 60)
        o = owner(stack)
        name = o === nothing ? NO_OWNER : o[2]
        line = o === nothing ? NO_OWNER : string(basename(o[1].file), ":", o[1].line, "  ", name)
        add!(total, type, a.size)
        any(fr -> fr.owned && fr.func == "time_loop!", stack) && add!(in_loop, type, a.size)
        add!(get!(AllocRow, by_func, name), type, a.size)
        add!(get!(AllocRow, by_line, line), type, a.size)
    end
    amount(row) = string(row.count, ", ", format_bytes(row.bytes))
    largest_first(rows) = sort!(collect(rows); by=kv -> (-last(kv).bytes, first(kv)))
    pct(row) = @sprintf("%.2f", 100*row.bytes/total.bytes)

    report, header = open_report("allocs", description, out_dir)
    push!(header, info..., "allocations" => amount(total),
        "in loop" => amount(in_loop) * " (with `time_loop!` on the stack)")

    close_report(report, "Allocations", header) do io, top
        rows_of(items) = top === nothing ? items : first(items, top)
        println(io, "## By function\n")
        println(io, "Every allocation of one call, charged to the innermost package frame of its stack.\n")
        md_table(io, ["Function", "Bytes", "Share %", "Count"],
            ([code(name), format_bytes(row.bytes), pct(row), string(row.count)] for (name, row) in rows_of(largest_first(by_func))))
        println(io, "## By line\n")
        println(io, "The ", top === nothing ? 40 : top, " largest, with the type holding the most bytes there.\n")
        md_table(io, ["Line", "Bytes", "Share %", "Count", "Largest type"],
            ([code(line), format_bytes(row.bytes), pct(row), string(row.count), first(largest(row.types, 1), 80)]
            for (line, row) in first(largest_first(by_line), top === nothing ? 40 : top)); left=(1, 5))
    end
    return report
end

"""
Check `f` for runtime dispatch with JET (needs JET.jl in the default environment;
static, nothing is timed); returns the report directory. Only calls in the package are reported.
"""
function dispatch(f; description::AbstractString="a closure", info=Pair{String,String}[], out_dir=nothing)
    JET = load_optional("JET", "c3a54625-cd67-489e-a8e7-0a5a0ff4e31b")   # before the analysis, so a missing package fails early
    result = Base.invokelatest(JET.report_opt, f, (); target_modules=(GREBClimate,))
    n = length(Base.invokelatest(JET.get_reports, result))

    report, header = open_report("dispatch", description, out_dir)
    push!(header, info..., "reports" => string(n))
    close_report(report, "Runtime dispatch", header) do io, _
        println(io, "JET optimization analysis of ", description, ", reports from the package only. ",
            "A report is a call that is dispatched at run time; none is the expected state.\n")
        if n > 0
            println(io, "```")
            Base.invokelatest(show, IOContext(io, :color => false), MIME"text/plain"(), result)
            println(io, "\n```")
        end
    end
    return report
end

# ----- Comparing two runs -----

"A run's report directory: `path` itself, or the folder of that name under `benchmark/profiles/`."
function run_dir(path::AbstractString)
    isdir(path) && return path
    inside = joinpath(PROFILES_DIR, path)
    isdir(inside) || error("no report directory at $path or $inside")
    return inside
end

"Read a run's `counts.tsv` as `(commit, samples, functions, lines, classes)`."
function read_run(dir::AbstractString)
    path = joinpath(dir, "counts.tsv")
    isfile(path) || error("$dir has no counts.tsv; run that profile again")
    tables = Dict(name => Dict{String,String}() for name in ("run", "class", "function", "line"))
    for row in Iterators.drop(eachline(path), 1)
        name, key, value = split(row, '\t')
        tables[name][key] = value
    end
    counts(name) = Dict(k => parse(Int, v) for (k, v) in tables[name])
    return (; commit=tables["run"]["commit"], samples=parse(Int, tables["run"]["samples"]),
        functions=counts("function"), lines=counts("line"), classes=counts("class"))
end

"Write one table of `compare_runs`: each row's share in runs `a` and `b`, the difference and its standard error."
function write_changes(io::IO, title::AbstractString, a::Dict{String,Int}, na::Int, b::Dict{String,Int}, nb::Int;
        top::Int=25)
    rows = map(collect(union(keys(a), keys(b)))) do name
        pa, pb = get(a, name, 0) / na, get(b, name, 0) / nb
        se = sqrt(pa * (1 - pa) / na + pb * (1 - pb) / nb)
        (; name, pa, pb, diff=pb - pa, se)
    end
    sort!(rows; by=r -> (-abs(r.diff), r.name))
    println(io, "## ", title, "\n")
    md_table(io, ["Row", "A %", "B %", "Diff", "+/-", "3 SE"],
        ([code(r.name), @sprintf("%.2f", 100r.pa), @sprintf("%.2f", 100r.pb), @sprintf("%+.2f", 100r.diff),
             @sprintf("%.2f", 100r.se), abs(r.diff) >= 3 * r.se ? "*" : ""] for r in first(rows, top)))
end

"Compare the shares of two runs, largest change first."
function compare_runs(dir_a::AbstractString, dir_b::AbstractString; io::IO=stdout)
    a, b = read_run(run_dir(dir_a)), read_run(run_dir(dir_b))
    println(io, "# Profile comparison\n")
    md_table(io, ["Run", "Directory", "Commit", "Samples"],
        ([["A", dir_a, a.commit, string(a.samples)], ["B", dir_b, b.commit, string(b.samples)]]); left=(1, 2, 3))
    println(io, "Self shares in percent; Diff is B minus A in percentage points, +/- its standard error. ",
        "* marks a difference of three standard errors or more. A row one run lacks counts as zero there.\n")
    write_changes(io, "By class", a.classes, a.samples, b.classes, b.samples)
    write_changes(io, "By function", a.functions, a.samples, b.functions, b.samples)
    write_changes(io, "By line", a.lines, a.samples, b.lines, b.samples)
    return nothing
end

end # module
