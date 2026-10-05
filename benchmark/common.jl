# Helpers shared by the benchmark scripts.

const REPO = normpath(joinpath(@__DIR__, ".."))
const JULIA_BIN = joinpath(Sys.BINDIR, Base.julia_exename())

"The local dataset directory; never downloads."
function default_data_dir()
    resolved = try
        greb_data_dir(; allow_download=false)
    catch err
        @warn "greb_data_dir could not resolve a dataset; falling back to the repo-local path" err
        nothing
    end
    something(resolved, joinpath(REPO, "greb_input_data"))
end

"Parse a repetition count (an integer >= 1)."
function parse_reps(s::AbstractString)
    n = tryparse(Int, s)
    n === nothing && error("reps must be an integer, got $(repr(s))")
    n >= 1 || error("reps must be at least 1, got $n")
    n
end

"Parse a non-negative integer flag value, e.g. for `--ctrl=10`."
function parse_nonneg_int(name::AbstractString, s::AbstractString)
    n = tryparse(Int, s)
    n === nothing && error("$name must be an integer, got $(repr(s))")
    n >= 0 || error("$name must be >= 0, got $n")
    n
end

"""
Split `--name=value` flags out of `args`. Returns `(flags, positional)`, a
`Dict{String,String}` of flag name to value and the remaining args in order.
Flags can appear anywhere; they don't have to come first.
"""
function split_flags(args::Vector{String})
    flags = Dict{String,String}()
    positional = String[]
    for a in args
        m = match(r"^--([\w-]+)=(.*)$", a)
        if m === nothing
            push!(positional, a)
        else
            flags[m.captures[1]] = m.captures[2]
        end
    end
    return flags, positional
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
    return "load average $(round(Sys.loadavg()[1], digits=2)) on $(Sys.CPU_THREADS) threads"
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
the calibration loop against the best seen and against `before`, and for
`times` the ratio of the median to the minimum.
"""
function machine_footer(before::Float64, times::Vector{Float64}=Float64[])
    after = calibrate()
    best = record_calibration(after)
    worst = max(before, after) / best
    reasons = String[]
    worst > NOISE_LIMIT && push!(reasons, "calibration $(round(worst, digits=2))x the best")
    if length(times) >= 3
        sorted = sort(times)
        spread = sorted[(end + 1) ÷ 2] / sorted[1]
        println("median / minimum of the runs: ", round(spread, digits=2))
        spread > NOISE_LIMIT && push!(reasons, "runs spread $(round(spread, digits=2))x")
    end
    println("calibration loop after: ", round(1e3 * after, digits=2), " ms, ", round(after / best, digits=2), "x the best")
    println(isempty(reasons) ? "machine state: quiet" : "machine state: NOISY (" * join(reasons, "; ") * ") - do not quote these timings")
    return isempty(reasons)
end
