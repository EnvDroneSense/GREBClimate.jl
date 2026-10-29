# Timing and allocation benchmarks for GREBClimate.jl.
#
#   julia --project=benchmark -t 2,0 benchmark/run_benchmarks.jl [mode] [jld2_dir] [reps]
#
#   year    - time a 1-year control run (default); --out=PATH writes the run times to a file
#   stages  - time each physics stage of one timestep; --seconds=S per stage (default 0.5)
#   threads - time `year` at -t 1, 2, 3, 4
#   alloc   - bytes allocated by one tendencies! call
#   years   - time a multi-year control+scenario run, e.g.:
#               julia --project=benchmark -t 2,0 benchmark/run_benchmarks.jl years --ctrl=10 --scnr=100
#             --ctrl=N, --scnr=N (default 10/10), --experiment=NAME (default full_model)
#   compare - circulation!(Ta) against circulation!(q), interleaved, with a self-check;
#             --seconds=S per variant (default 5); --counters adds instructions and
#             cycles per call (Linux, needs LinuxPerf in your own Julia environment). Exits non-zero when the self-check fails.
#             For other variants call compare_kernels (benchmark/compare.jl) from a script.
#
# Run it as `julia --project=benchmark` (Chairmarks), with the package added by
# `Pkg.develop(path=".")` in that environment.

using GREBClimate
using GREBClimate: SWradiation!, LWradiation!, hydro!, convergence!, seaice!, deep_ocean!,
                   circulation!, tendencies!, output!, diagnostics!, init_model!

include("common.jl")
include("compare.jl")

"True if the dataset exists; warns otherwise (benchmarks never download it)."
function _require_data(jld2_dir::AbstractString)
    isdir(jld2_dir) && return true
    @warn "JLD2 data directory not found: $jld2_dir. Set GREB_DATA or pass a path."
    return false
end

"Run `greb_model!` without its log output."
quiet_run(spec, cfg, jld2_dir, fields) = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
    greb_model!(spec, cfg; jld2_dir=jld2_dir, fields=fields)
end

"Seconds for one `quiet_run` on a fresh copy of `fields`."
time_run(spec, cfg, jld2_dir, fields) = (@b deepcopy(fields) quiet_run(spec, cfg, jld2_dir, _) evals=1 samples=1).time

"Time a 1-year `:full_model` control run (stored corrections, no spin-up) `reps` times; returns seconds per run."
function time_1yr(jld2_dir::AbstractString; cfg=preset(:full_model; corrections=Stored()), reps::Int=3)
    _require_data(jld2_dir) || return nothing
    reps >= 1 || throw(ArgumentError("reps must be at least 1, got $reps"))

    println("Threads.nthreads() = ", Threads.nthreads())
    calibration = machine_header()
    fields = load_climatology(jld2_dir; dataset=:ncep)

    # Warm-up, so compilation is not timed.
    quiet_run(RunSpec(scnr=0), cfg, jld2_dir, deepcopy(fields))

    times = Float64[]
    for r in 1:reps
        push!(times, time_run(RunSpec(scnr=0), cfg, jld2_dir, fields))
        println("  run $r: ", round(times[end], digits=3), " s")
    end

    println("min: ", round(minimum(times), digits=3), " s  median: ", round(median(times), digits=3), " s")
    record_result("year", "year", times, machine_footer(calibration, times))
    return times
end

"""
Time a `ctrl`-year control run followed by a `scnr`-year scenario run, `reps`
times on the stored flux corrections (no spin-up); returns seconds per run.
Unlike `time_1yr` (fixed 1-year control run, no scenario), this is for checking long-run cost and stability.

`experiment` is a preset name. `:full_model` holds CO2 constant across the
whole scenario regardless of `scnr`; for an actual multi-year forced change,
pass a preset with a real CO2 trajectory.
"""
function time_years(jld2_dir::AbstractString; experiment::Symbol=:full_model,
        ctrl::Int=10, scnr::Int=10, reps::Int=1)
    _require_data(jld2_dir) || return nothing
    reps >= 1 || throw(ArgumentError("reps must be at least 1, got $reps"))
    ctrl >= 1 || throw(ArgumentError("ctrl must be at least 1, got $ctrl"))
    scnr >= 0 || throw(ArgumentError("scnr must be >= 0, got $scnr"))
    total_years = ctrl + scnr

    println("Threads.nthreads() = ", Threads.nthreads())
    println("ctrl=$ctrl scnr=$scnr ($total_years simulated years/rep), experiment=$experiment")
    calibration = machine_header()
    cfg = preset(experiment; corrections=Stored())  # no spin-up, as in time_1yr
    fields = load_climatology(jld2_dir; dataset=:ncep)

    # Warm-up with a minimal run.
    quiet_run(RunSpec(ctrl=1, scnr=0), cfg, jld2_dir, deepcopy(fields))

    times = Float64[]
    for r in 1:reps
        push!(times, time_run(RunSpec(ctrl=ctrl, scnr=scnr), cfg, jld2_dir, fields))
        println("  run $r: ", round(times[end], digits=3), " s  (",
            round(times[end] / total_years, digits=3), " s/simulated year)")
    end

    println("min: ", round(minimum(times), digits=3), " s  median: ", round(median(times), digits=3), " s  ",
        "= ", round(median(times) / total_years, digits=3), " s/simulated year (median) over ", total_years, " years")
    record_result("years", experiment, times, machine_footer(calibration, times))
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

"""
Time each physics stage of one timestep for about `seconds` each; returns
`(name, fastest seconds per call)`. The stage is timed in `ROUNDS` chunks, and
the median is taken over the fastest call of each chunk.
"""
function time_stages(jld2_dir::AbstractString; cfg=preset(:full_model), seconds::Real=0.5)
    _require_data(jld2_dir) || return nothing
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
        println("  ", rpad(name, 18), "min ", round(per_call * 1e6, digits=2), " µs   median ",
            round(median(t) * 1e6, digits=2), " µs   (", round(share, digits=1), "% of measured total)")
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
    _require_data(jld2_dir) || return nothing

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
    if haskey(results, base_n)
        base = sum(results[base_n]) / length(results[base_n])
        println("\nSpeedup vs -t $base_n:")
        for n in thread_counts
            haskey(results, n) || continue
            r = results[n]
            mean_n = sum(r) / length(r)
            println("  -t $n: ", round(mean_n, digits=3), "s  ",
                "(min ", round(minimum(r), digits=3), "s, max ", round(maximum(r), digits=3), "s)  ",
                "(", round(base / mean_n, digits=2), "x)")
        end
    end
    return results
end

"""
Compare `circulation!(Ta)` with `circulation!(q)` (the baseline is the first);
returns whether the self-check passed. A template for `compare_kernels`.
"""
function compare_circulation(jld2_dir::AbstractString; seconds::Real=5, counters::Bool=false)
    _require_data(jld2_dir) || return nothing
    (; fields, r, ws, slice, Ta, q) = stage_inputs(jld2_dir, preset(:full_model))
    p = r.config.processes
    return compare_kernels(
        "circulation!(Ta)" => () -> circulation!(ws.dTa_crcl, Ta, GREBClimate.z_air, fields, ws, slice, p),
        "circulation!(q)" => () -> circulation!(ws.dq_crcl, q, GREBClimate.z_vapor, fields, ws, slice, p);
        seconds, counters)
end

const TENDENCIES_ALLOC_BUDGET = 256  # kept in sync with test/test_invariants.jl

"Bytes allocated by one `tendencies!` call, checked against the test budget."
function check_allocations(jld2_dir::AbstractString)
    _require_data(jld2_dir) || return nothing

    (; fields, r, CO2, state, ws, slice, Ts, Ta, To, q) = stage_inputs(jld2_dir, preset(:full_model))

    tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, slice, r)  # warm-up
    bytes = @allocated tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, slice, r)

    verdict = bytes <= TENDENCIES_ALLOC_BUDGET ? "within" : "OVER"
    println("tendencies! allocations (single-workspace path): ", bytes, " bytes ",
        "($verdict the $TENDENCIES_ALLOC_BUDGET-byte budget in test/test_invariants.jl)")
    return bytes
end

const _MODES = ("year", "stages", "threads", "alloc", "years", "compare")

if abspath(PROGRAM_FILE) == @__FILE__
    mode, rest = if isempty(ARGS)
        ("year", String[])
    elseif ARGS[1] in _MODES
        (ARGS[1], ARGS[2:end])
    else
        error("unknown mode $(repr(ARGS[1])); expected one of $(join(_MODES, ", "))")
    end

    # `--name=value` flags can appear anywhere in `rest`, so they're stripped
    # before the positional [jld2_dir] [reps].
    flags, rest = split_flags(rest)
    counters = "--counters" in rest
    filter!(!=("--counters"), rest)

    jld2_dir = !isempty(rest) ? rest[1] : default_data_dir()

    reps = length(rest) >= 2 ? parse_reps(rest[2]) : nothing

    seconds = haskey(flags, "seconds") ? parse(Float64, flags["seconds"]) : nothing

    if mode == "year"
        times = time_1yr(jld2_dir; reps=something(reps, 3))
        haskey(flags, "out") && times !== nothing && write(flags["out"], join(times, "\n"))
    elseif mode == "stages"
        reps === nothing || error("the stages mode takes --seconds=S, not a reps argument")
        time_stages(jld2_dir; seconds=something(seconds, 0.5))
    elseif mode == "compare"
        reps === nothing || error("the compare mode takes --seconds=S, not a reps argument")
        compare_circulation(jld2_dir; seconds=something(seconds, 5), counters) === false && exit(1)
    elseif mode == "threads"
        sweep_threads(jld2_dir; reps=something(reps, 3))
    elseif mode == "alloc"
        reps === nothing || error("the alloc mode takes no reps argument")
        check_allocations(jld2_dir)
    elseif mode == "years"
        ctrl = haskey(flags, "ctrl") ? parse_nonneg_int("ctrl", flags["ctrl"]) : 10
        scnr = haskey(flags, "scnr") ? parse_nonneg_int("scnr", flags["scnr"]) : 10
        experiment = Symbol(get(flags, "experiment", "full_model"))
        time_years(jld2_dir; experiment,
            ctrl=ctrl, scnr=scnr, reps=something(reps, 1))
    end
end
