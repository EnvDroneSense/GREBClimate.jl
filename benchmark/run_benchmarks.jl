# =============================================================================
# run_benchmarks.jl - dependency-free timing harness for GREBClimate.jl
#
# Modes (see each function's docstring for what it measures):
#   year     - wall-clock for a 1-simulated-year control run  -> time_1yr
#   stages   - per-timestep breakdown by kernel                -> time_stages
#   threads  - `year` across -t N subprocesses                 -> sweep_threads
#   alloc    - bytes allocated by one tendencies! call         -> check_allocations
#
#   julia --project=. -t 2 benchmark/run_benchmarks.jl [mode] [jld2_dir] [reps]
#
# mode defaults to "year"; jld2_dir is resolved by `greb_data_dir` ($GREB_DATA
# or ../greb_input_data, never a download); `reps` applies to year/stages/
# threads, each with its own default. `-t 2`, not `-t 3`, is the recommended
# count - see the `sweep_threads` docstring. Also includable from a REPL:
# `include("benchmark/run_benchmarks.jl"); time_1yr(dir)`.
#
# `fields` and `cfg` are deep-copied per repetition, always OUTSIDE the timed
# region: `greb_model!` mutates `fields` in place and `init_model!` writes the
# hydrology parameters and any scenario table back onto `cfg`. Copying is a
# large thread-count-independent constant that would otherwise swamp the
# signal. `.claude/skills/benchmark/SKILL.md` has this machine's noise sources.
# =============================================================================

using GREBClimate

function _require_data(jld2_dir::AbstractString)
    isdir(jld2_dir) && return true
    @warn """
          JLD2 data directory not found: $jld2_dir

          Benchmarks never download data (resolved with allow_download=false).
          Fetch it once via `julia --project=. examples/run_greb.jl`, or point
          at an existing copy with GREB_DATA or a positional argument.
          """
    return false
end

"""
    time_1yr(jld2_dir; cfg=create_experiment_config(:full_model), reps=3)

Runs a real 1-simulated-year `:full_model` control run `reps` times (after
one untimed warm-up call to exclude JIT compilation) and prints each timing
plus the mean/min/max. Returns the vector of timings in seconds, or `nothing`
if the dataset is missing.
"""
function time_1yr(jld2_dir::AbstractString; cfg=create_experiment_config(:full_model), reps::Int=3)
    _require_data(jld2_dir) || return nothing
    reps >= 1 || throw(ArgumentError("reps must be at least 1, got $reps"))

    println("Threads.nthreads() = ", Threads.nthreads())
    fields = load_greb_jld2!(jld2_dir; dataset=:ncep)

    # Warm-up: excludes JIT compilation from every timed run below.
    redirect_stdout(devnull) do
        greb_model!(RunSpec(scnr=0), deepcopy(cfg); jld2_dir=jld2_dir, fields=deepcopy(fields))
    end

    times = Float64[]
    for r in 1:reps
        fields_r = deepcopy(fields)  # copies stay outside the timed region
        cfg_r = deepcopy(cfg)
        t = @elapsed redirect_stdout(devnull) do
            greb_model!(RunSpec(scnr=0), cfg_r; jld2_dir=jld2_dir, fields=fields_r)
        end
        push!(times, t)
        println("  run $r: ", round(t, digits=3), " s")
    end

    println("mean: ", round(sum(times) / length(times), digits=3), " s  ",
        "(min ", round(minimum(times), digits=3), "s, max ", round(maximum(times), digits=3), "s)")
    return times
end

"""
    time_stages(jld2_dir; cfg=create_experiment_config(:full_model), reps=2000)

Breaks one timestep's `tendencies!` pipeline into its component stages
(circulation for `Ta` and `q`, SW/LW radiation, `hydro!`, `deep_ocean!`) and
times each individually - `reps` back-to-back calls after a warm-up,
averaged. Returns `(name, seconds_per_call)` pairs, or `nothing` if the
dataset is missing.

`circulation!` is not a leaf: it runs `diffusion!`/`advection!` `ntime` times,
plus `convergence!` for `q` only (the moisture scale height selects it). So
`convergence!` is *included* here, and it is why the `q` lane costs more.

Single-workspace calls only, so nothing is spawned even at `-t 2`: this is the
serial cost that `tendencies!`'s parallel path races against.
"""
function time_stages(jld2_dir::AbstractString; cfg=create_experiment_config(:full_model), reps::Int=2000)
    _require_data(jld2_dir) || return nothing
    reps >= 1 || throw(ArgumentError("reps must be at least 1, got $reps"))

    fields = load_greb_jld2!(jld2_dir; dataset=:ncep)
    init_model!(cfg, fields)
    state = ModelState()
    ws = CirculationWorkspace()
    timestate = TimeState(1, 1)
    CO2 = cfg.co2_concentration

    ityr = timestate.ityr
    Ts = copy(fields.Tclim[:, :, ityr])
    Ta = copy(Ts)
    To = copy(fields.Toclim[:, :, ityr])
    q = copy(fields.qclim[:, :, ityr])

    stages = [
        ("circulation!(Ta)", () -> circulation!(Ta, GREBClimate.z_air, ws.dTa_crcl, fields, ws, timestate, cfg)),
        ("circulation!(q)", () -> circulation!(q, GREBClimate.z_vapor, ws.dq_crcl, fields, ws, timestate, cfg)),
        ("SWradiation!", () -> SWradiation!(Ts, fields, state, timestate, cfg, ws)),
        ("LWradiation!", () -> LWradiation!(Ts, Ta, q, CO2, fields, timestate, cfg, ws)),
        ("hydro!", () -> hydro!(Ts, q, fields, timestate, cfg, ws)),
        ("deep_ocean!", () -> deep_ocean!(Ts, To, fields, timestate, cfg, ws)),
    ]

    for (_, f) in stages
        f()
    end

    println("Per-stage timing (", reps, " calls each, single workspace, ",
        Threads.nthreads(), " thread(s) available but unused here):")
    results = Tuple{String,Float64}[]
    for (name, f) in stages
        t = @elapsed for _ in 1:reps
            f()
        end
        push!(results, (name, t / reps))
    end

    total = sum(last(r) for r in results)
    for (name, per_call) in results
        share = 100 * per_call / total
        println("  ", rpad(name, 18), round(per_call * 1e6, digits=2), " µs/call   (",
            round(share, digits=1), "% of measured total)")
    end
    println("  measured total (sum of stages): ", round(total * 1e6, digits=1), " µs")
    println("  (convergence! is inside circulation!(q); seaice!/output!/diagnostics!")
    println("   and the tendency assembly are not measured - see time_1yr for whole-model cost)")
    return results
end

"""
    sweep_threads(jld2_dir; thread_counts=(1,2,3,4), reps=3)

Runs `time_1yr` in separate `julia -t N` subprocesses for each `N` in
`thread_counts` - thread count is fixed at Julia startup, so this can't be
swept in-process - and reports each run's mean and range plus relative
speedup vs. the first entry. Returns a `Dict` mapping thread count to that
subprocess's vector of per-run timings.

Why `-t 2` and not more: with `nthreads() > 1` and distinct workspaces,
`tendencies!` `@spawn`s `circulation!(Ta)` and `circulation!(q)` and keeps the
rest on the calling task - but circulation is 98.7% of the stage total (see
`time_stages`), leaving ~30µs of synchronous work that cannot fill a third
thread. Measured over 5 sweeps: `-t 2` a consistent ~1.5× (1.24-1.62×);
`-t 3`/`-t 4` a noisy 1.0-1.65×, no dependable edge. Only `-t 2` vs `-t 1`
holds up in a single run.

Two limits. This is the *control* phase only - `RunSpec(scnr=0)` runs no
flux-correction spin-up, and `greb_model!` calls `qflux_correction!` with
`ws_a = ws_q = ws` (src/model.jl), pinning that phase to the serial path at
any `-t N`, so a flux-corrected run scales worse than these numbers. And
every subprocess reloads the climatology, which dominates the sweep's
wall-clock though not the timings.
"""
function sweep_threads(jld2_dir::AbstractString; thread_counts=(1, 2, 3, 4), reps::Int=3)
    _require_data(jld2_dir) || return nothing

    project_dir = joinpath(@__DIR__, "..")
    julia_bin = joinpath(Sys.BINDIR, Base.julia_exename())
    script = @__FILE__

    results = Dict{Int,Vector{Float64}}()
    for n in thread_counts
        println("--- -t $n ---")
        cmd = `$julia_bin --project=$project_dir -t $n $script year $jld2_dir $reps`
        output = read(cmd, String)
        print(output)
        # Per-run timings, not just the mean: a single number hides the spread.
        runs = [parse(Float64, m.captures[1]) for m in eachmatch(r"run\s+\d+:\s*([\d.]+)\s*s", output)]
        if isempty(runs)
            @warn "Could not parse any run timings from -t $n run"
        else
            results[n] = runs
        end
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

const TENDENCIES_ALLOC_BUDGET = 256  # kept in sync with test/test_invariants.jl

"""
    check_allocations(jld2_dir)

Reports bytes allocated by one `tendencies!` call (after a warm-up call to
exclude JIT/compilation allocations) and whether it is within
`TENDENCIES_ALLOC_BUDGET`. Uses the default `ws_a=ws_q=ws` (synchronous) path,
since spawning tasks allocates regardless of the physics code.

Reports only; `test/test_invariants.jl` is what *enforces* the budget, for
every kernel, on every CI run.
"""
function check_allocations(jld2_dir::AbstractString)
    _require_data(jld2_dir) || return nothing

    cfg = create_experiment_config(:full_model)
    fields = load_greb_jld2!(jld2_dir; dataset=:ncep)
    init_model!(cfg, fields)
    state = ModelState()
    ws = CirculationWorkspace()
    timestate = TimeState(1, 1)
    CO2 = cfg.co2_concentration

    ityr = timestate.ityr
    Ts = copy(fields.Tclim[:, :, ityr])
    Ta = copy(Ts)
    To = copy(fields.Toclim[:, :, ityr])
    q = copy(fields.qclim[:, :, ityr])

    tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, timestate, cfg)  # warm-up
    bytes = @allocated tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, timestate, cfg)

    verdict = bytes <= TENDENCIES_ALLOC_BUDGET ? "within" : "OVER"
    println("tendencies! allocations (single-workspace path): ", bytes, " bytes ",
        "($verdict the $TENDENCIES_ALLOC_BUDGET-byte budget in test/test_invariants.jl)")
    return bytes
end

"""
    default_data_dir() -> String

Dataset directory to use when none was given on the command line.

A function, not a top-level `const`: `greb_data_dir` raises on a `GREB_DATA`
that no longer exists, and evaluating that at load time broke every mode -
even ones handed an explicit path, even `include` from a REPL.
`allow_download=false` so benchmarking can never pull 353 MB as a side effect.
"""
function default_data_dir()
    resolved = try
        greb_data_dir(; allow_download=false)
    catch err
        @warn "greb_data_dir could not resolve a dataset; falling back to the repo-local path" err
        nothing
    end
    return something(resolved, joinpath(@__DIR__, "..", "greb_input_data"))
end

const _MODES = ("year", "stages", "threads", "alloc")

# Runs as a script but not when `include`-d, so a REPL is never terminated.
# ARGS[1] is taken as a directory (the legacy call form) only when it names
# one, so a mistyped mode reports itself instead of a missing dataset.
if abspath(PROGRAM_FILE) == @__FILE__
    mode, rest = if isempty(ARGS)
        ("year", String[])
    elseif ARGS[1] in _MODES
        (ARGS[1], ARGS[2:end])
    elseif isdir(ARGS[1])
        ("year", ARGS)          # legacy form: bare directory, no mode
    else
        error("unknown mode $(repr(ARGS[1])); expected one of $(join(_MODES, ", ")) " *
              "or a path to an existing dataset directory")
    end

    jld2_dir = !isempty(rest) ? rest[1] : default_data_dir()

    reps = if length(rest) >= 2
        parsed = tryparse(Int, rest[2])
        parsed === nothing && error("reps must be an integer, got $(repr(rest[2]))")
        parsed < 1 && error("reps must be at least 1, got $parsed")
        parsed
    else
        nothing
    end

    if mode == "year"
        time_1yr(jld2_dir; reps=something(reps, 3))
    elseif mode == "stages"
        time_stages(jld2_dir; reps=something(reps, 2000))
    elseif mode == "threads"
        sweep_threads(jld2_dir; reps=something(reps, 3))
    elseif mode == "alloc"
        # One measurement; accepting reps would silently do nothing.
        reps === nothing || error("the alloc mode takes no reps argument")
        check_allocations(jld2_dir)
    end
end
