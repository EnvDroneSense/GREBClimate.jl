# Timing benchmarks for GREBClimate.jl.
#
#   julia --project=benchmark -t 2,0 benchmark/run_benchmarks.jl [mode] [jld2_dir] [reps]
#
#   year    - time a 1-year control run (default); --out=PATH writes the run times to a file
#   stages  - time each physics stage of one timestep; --seconds=S per stage (default 0.5)
#   threads - time `year` at -t 1, 2, 3, 4
#   years   - time a multi-year control+scenario run, e.g.:
#               julia --project=benchmark -t 2,0 benchmark/run_benchmarks.jl years --ctrl=10 --scnr=100
#             --ctrl=N, --scnr=N (default 10/10), --experiment=NAME (default full_model)
#   compare - circulation!(Ta) against circulation!(q), interleaved, with a self-check;
#             --seconds=S per variant (default 5); --counters adds instructions and
#             cycles per call (Linux, needs LinuxPerf in your own Julia environment). Exits non-zero when the self-check fails.
#             For other variants call compare_kernels (benchmark/compare.jl) from a script.
#
# A flag a mode does not take, or a missing dataset, stops the run with exit code 1.
# Run it as `julia --project=benchmark` (Chairmarks), with the package added by
# `Pkg.develop(path=".")` in that environment. Allocations are checked by
# test/test_invariants.jl, not here.

using GREBClimate
using GREBClimate: SWradiation!, LWradiation!, hydro!, convergence!, seaice!, deep_ocean!,
                   circulation!, tendencies!, output!, diagnostics!, init_model!

include("common.jl")
include("compare.jl")

"Run `greb_model!` without its log output."
quiet_run(spec, cfg, jld2_dir, fields) = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
    greb_model!(spec, cfg; jld2_dir=jld2_dir, fields=fields)
end

"""
Time a `ctrl`-year control run followed by a `scnr`-year scenario run, `reps`
times on the stored flux corrections (no spin-up); returns seconds per run.
`mode` is the name in the result record: `"year"` for the 1-year control run
of the `year` mode, `"years"` otherwise.

`experiment` is a preset name. `:full_model` holds CO2 constant across the
whole scenario regardless of `scnr`; for an actual multi-year forced change,
pass a preset with a real CO2 trajectory.
"""
function time_years(jld2_dir::AbstractString; mode::AbstractString="years", experiment::Symbol=:full_model,
        ctrl::Int=10, scnr::Int=10, reps::Int=1)
    reps >= 1 || throw(ArgumentError("reps must be at least 1, got $reps"))
    ctrl >= 1 || throw(ArgumentError("ctrl must be at least 1, got $ctrl"))
    scnr >= 0 || throw(ArgumentError("scnr must be >= 0, got $scnr"))
    total_years = ctrl + scnr

    println("Threads.nthreads() = ", Threads.nthreads())
    println("ctrl=$ctrl scnr=$scnr ($total_years simulated years/rep), experiment=$experiment")
    calibration = machine_header()
    cfg = preset(experiment; corrections=Stored())  # stored corrections: no spin-up
    fields = load_climatology(jld2_dir; dataset=:ncep)
    spec = RunSpec(ctrl=ctrl, scnr=scnr)

    # Warm-up with a minimal run, so compilation is not timed. `greb_model!`
    # restores what it mutates in `fields`, so the same copy serves every run.
    quiet_run(RunSpec(ctrl=1, scnr=0), cfg, jld2_dir, fields)

    times = Float64[]
    for r in 1:reps
        push!(times, @elapsed quiet_run(spec, cfg, jld2_dir, fields))
        println("  run $r: ", round(times[end], digits=3), " s  (",
            round(times[end] / total_years, digits=3), " s/simulated year)")
    end

    println("min: ", round(minimum(times), digits=3), " s  median: ", round(_median(times), digits=3), " s  ",
        "= ", round(_median(times) / total_years, digits=3), " s/simulated year (median) over ", total_years, " years")
    record_result(mode, mode == "year" ? "year" : String(experiment), times, machine_footer(calibration, times))
    return times
end

"The climatology, workspace and slice-1 state that the single-stage modes call the kernels with."
function stage_inputs(jld2_dir::AbstractString, cfg)
    fields = load_climatology(jld2_dir; dataset=:ncep)
    r = resolve(cfg; jld2_dir)
    CO2 = init_model!(fields, r).CO2_ctrl
    slice = 1
    Ts = copy(fields.Ts_clim[:, :, slice])
    return (; fields, r, CO2, state=ModelState(), ws=ModelWorkspace(), slice, Ts, Ta=copy(Ts),
        To=copy(fields.To_clim[:, :, slice]), q=copy(fields.q_clim[:, :, slice]))
end

"`x` with `digits` decimals (trailing zeros kept), right-aligned in `width` characters."
function _fixed(x::Real, digits::Int, width::Int)
    whole, frac = split(string(round(x; digits)), '.')
    return lpad(whole * "." * rpad(frac, digits, '0'), width)
end

"""
Time each physics stage of one timestep for about `seconds` each; returns
`(name, fastest seconds per call)`. The stage is timed in `ROUNDS` chunks; the
fastest call of each chunk is kept, and the minimum and median of those are printed.
"""
function time_stages(jld2_dir::AbstractString; cfg=preset(:full_model), seconds::Real=0.5)
    seconds > 0 || throw(ArgumentError("seconds must be positive, got $seconds"))

    (; fields, r, CO2, state, ws, slice, Ts, Ta, To, q) = stage_inputs(jld2_dir, cfg)
    p = r.config.processes

    stages = [
        ("circulation!(Ta)", () -> circulation!(ws.dTa_crcl, Ta, GREBClimate.z_air, fields, ws, slice, p)),
        ("circulation!(q)", () -> circulation!(ws.dq_crcl, q, GREBClimate.z_vapor, fields, ws, slice, p)),
        ("SWradiation!", () -> SWradiation!(ws, Ts, fields, state, slice, p)),
        ("LWradiation!", () -> LWradiation!(ws, Ts, Ta, q, CO2, fields, slice, p)),
        ("hydro!", () -> hydro!(ws, Ts, q, fields, slice, p, r.hydrology)),
        ("deep_ocean!", () -> deep_ocean!(ws, Ts, To, fields, slice, p))
    ]

    for (_, f) in stages
        f()
    end
    calibration = machine_header()

    println("Per-stage timing (about ", seconds, " s each, single workspace, ",
        Threads.nthreads(), " thread(s) available but unused here):")
    samples = [[chunk_time(f, seconds / ROUNDS) for _ in 1:ROUNDS] for (_, f) in stages]
    results = [(name, minimum(t)) for ((name, _), t) in zip(stages, samples)]

    total = sum(last(r) for r in results)
    for ((name, per_call), t) in zip(results, samples)
        share = 100 * per_call / total
        println("  ", rpad(name, 18), "min ", _fixed(per_call * 1e6, 2, 7), " µs   median ",
            _fixed(_median(t) * 1e6, 2, 7), " µs   (", _fixed(share, 1, 5), "% of measured total)")
    end
    println("  measured total (sum of stage minima): ", round(total * 1e6, digits=1), " µs")
    println("  (convergence! is inside circulation!(q). The shares are of these six calls, not of")
    println("   a step: seaice!, the update loop, output! and diagnostics! are not timed. For shares")
    println("   of a whole run use benchmark/profile.jl step)")
    quiet = machine_footer(calibration)
    for ((name, _), t) in zip(stages, samples)
        record_result("stages", name, t, quiet)
    end
    return results
end

"Time `year` in a subprocess per thread count; returns seconds per run for each."
function sweep_threads(jld2_dir::AbstractString; thread_counts=(1, 2, 3, 4), reps::Int=3)
    script = @__FILE__
    project = dirname(Base.active_project())

    results = Dict{Int,Vector{Float64}}()
    for n in thread_counts
        println("--- -t $n ---")
        out = tempname()
        run(`$JULIA_BIN --project=$project -t $n,0 $script year $jld2_dir $reps --out=$out`)
        results[n] = parse.(Float64, readlines(out))
        rm(out)
    end

    base_n = first(thread_counts)
    base = _median(results[base_n])
    println("\nSpeedup of the median vs -t $base_n:")
    for n in thread_counts
        r = results[n]
        println("  -t $n: min ", round(minimum(r), digits=3), " s  median ", round(_median(r), digits=3),
            " s  (", round(base / _median(r), digits=2), "x)")
    end
    return results
end

"""
Compare `circulation!(Ta)` with `circulation!(q)` (the baseline is the first);
returns whether the self-check passed. A template for `compare_kernels`.
"""
function compare_circulation(jld2_dir::AbstractString; seconds::Real=5, counters::Bool=false)
    (; fields, r, ws, slice, Ta, q) = stage_inputs(jld2_dir, preset(:full_model))
    p = r.config.processes
    return compare_kernels(
        "circulation!(Ta)" => () -> circulation!(ws.dTa_crcl, Ta, GREBClimate.z_air, fields, ws, slice, p),
        "circulation!(q)" => () -> circulation!(ws.dq_crcl, q, GREBClimate.z_vapor, fields, ws, slice, p);
        seconds, counters)
end

const _MODES = ("year", "stages", "threads", "years", "compare")

# The flags each mode takes, and whether it takes the positional [reps].
const _MODE_FLAGS = (year=("out",), stages=("seconds",), threads=(), years=("ctrl", "scnr", "experiment"),
    compare=("seconds", "counters"))
const _MODE_REPS = (year=true, stages=false, threads=true, years=true, compare=false)

"Run the benchmark named by `args`; returns the exit code."
function main(args::Vector{String})
    mode, rest = if isempty(args)
        ("year", String[])
    elseif args[1] in _MODES
        (args[1], args[2:end])
    else
        error("unknown mode $(repr(args[1])); expected one of $(join(_MODES, ", "))")
    end
    sym = Symbol(mode)

    # Flags can appear anywhere in `rest`, so they are stripped before the
    # positional [jld2_dir] [reps].
    flags, rest = split_flags(rest)
    for name in keys(flags)
        name in _MODE_FLAGS[sym] || error("the $mode mode takes no --$name flag")
    end
    length(rest) <= (_MODE_REPS[sym] ? 2 : 1) ||
        error("the $mode mode takes ", _MODE_REPS[sym] ? "[jld2_dir] [reps]" : "[jld2_dir] and no reps argument", ", got: ", join(rest, " "))

    jld2_dir = !isempty(rest) ? rest[1] : default_data_dir()
    isdir(jld2_dir) || error("dataset directory not found: $jld2_dir. Set GREB_DATA or pass a path.")
    reps = length(rest) >= 2 ? parse_reps(rest[2]) : nothing

    if mode == "year"
        times = time_years(jld2_dir; mode, ctrl=1, scnr=0, reps=something(reps, 3))
        haskey(flags, "out") && write(flags["out"], join(times, "\n"))
    elseif mode == "stages"
        seconds = haskey(flags, "seconds") ? parse(Float64, flags["seconds"]) : 0.5
        time_stages(jld2_dir; seconds)
    elseif mode == "compare"
        seconds = haskey(flags, "seconds") ? parse(Float64, flags["seconds"]) : 5
        compare_circulation(jld2_dir; seconds, counters=haskey(flags, "counters")) || return 1
    elseif mode == "threads"
        sweep_threads(jld2_dir; reps=something(reps, 3))
    else  # years
        ctrl = haskey(flags, "ctrl") ? parse_nonneg_int("ctrl", flags["ctrl"]; min=1) : 10
        scnr = haskey(flags, "scnr") ? parse_nonneg_int("scnr", flags["scnr"]) : 10
        experiment = Symbol(get(flags, "experiment", "full_model"))
        time_years(jld2_dir; mode, experiment, ctrl, scnr, reps=something(reps, 1))
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
