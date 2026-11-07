# Helpers shared by the benchmark scripts.

using Dates

const REPO = normpath(joinpath(@__DIR__, ".."))
const JULIA_BIN = joinpath(Sys.BINDIR, Base.julia_exename())

"The local dataset directory; never downloads. Stops when there is none, or when `GREB_DATA` points nowhere."
function default_data_dir()
    dir = greb_data_dir(; allow_download=false)
    dir === nothing && error("no dataset found. Set GREB_DATA or pass the dataset directory.")
    return dir
end

"Parse an integer of at least `min` (default 0), e.g. for `--ctrl=10`."
function parse_nonneg_int(name::AbstractString, s::AbstractString; min::Int=0)
    n = tryparse(Int, s)
    n === nothing && error("$name must be an integer, got $(repr(s))")
    n >= min || error("$name must be at least $min, got $n")
    n
end

"Parse a repetition count (an integer >= 1)."
parse_reps(s::AbstractString) = parse_nonneg_int("reps", s; min=1)

"""
Split flags out of `args`. Returns `(flags, positional)`: a `Dict{String,String}`
of flag name to value (`--name=value`; a bare `--name` gets `"true"`) and the
remaining args in order. Flags can appear anywhere; they don't have to come first.
"""
function split_flags(args::Vector{String})
    flags = Dict{String,String}()
    positional = String[]
    for a in args
        m = match(r"^--([\w-]+)(?:=(.*))?$", a)
        if m === nothing
            push!(positional, a)
        else
            flags[m.captures[1]] = something(m.captures[2], "true")
        end
    end
    return flags, positional
end

"The median of `v`: the mean of the two middle values for an even length."
function _median(v::AbstractVector{<:Real})
    sorted = sort(v)
    return (sorted[(end+1)÷2] + sorted[(end+2)÷2]) / 2
end

# ----- Machine state -----
#
# A timing only means something next to how fast the machine was at that
# moment. A fixed loop is timed before and after the benchmark and compared
# with the best time this machine has ever given for it.

const CALIBRATION_FILE = joinpath(REPO, "benchmark", "profiles", "calibration.txt")
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

"Short commit hash of the repository and whether tracked files are modified."
function repo_state()
    try
        commit = readchomp(`git -C $REPO rev-parse --short HEAD`)
        dirty = !isempty(readchomp(`git -C $REPO status --porcelain --untracked-files=no`))
        return commit, dirty
    catch
        return "unknown", false
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
        spread = _median(times) / minimum(times)
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
ratio or `nothing`; `quiet` is what `machine_footer` returned.
"""
function record_result(mode::AbstractString, name::AbstractString, times::Vector{Float64}, quiet::Bool;
        ratios=nothing, selfcheck=nothing)
    commit, dirty = repo_state()
    fields = (Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"), commit * (dirty ? "-dirty" : ""), mode, name,
        Threads.nthreads(), length(times), minimum(times), _median(times),
        ratios === nothing ? "" : ratios[1], ratios === nothing ? "" : ratios[2],
        selfcheck === nothing ? "" : selfcheck, quiet ? "quiet" : "NOISY")
    mkpath(RESULTS_DIR)
    file = joinpath(RESULTS_DIR, Dates.format(today(), "yyyy-mm-dd") * ".csv")
    isfile(file) || write(file, RECORD_COLUMNS * "\n")
    open(io -> println(io, join(fields, ',')), file, "a")
    return nothing
end
