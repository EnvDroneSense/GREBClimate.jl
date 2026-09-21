# Timing and allocation benchmarks for GREBClimate.jl.
#
#   julia --project=. -t 2 benchmark/run_benchmarks.jl [mode] [jld2_dir] [reps]
#
#   year    - time a 1-year control run (default)
#   stages  - time each physics stage of one timestep
#   threads - time `year` at -t 1, 2, 3, 4
#   alloc   - bytes allocated by one tendencies! call

using GREBClimate

include("common.jl")

"True if the dataset exists; warns otherwise (benchmarks never download it)."
function _require_data(jld2_dir::AbstractString)
    isdir(jld2_dir) && return true
    @warn "JLD2 data directory not found: $jld2_dir. Set GREB_DATA or pass a path."
    return false
end

"Time a 1-year `:full_model` control run `reps` times; returns seconds per run."
function time_1yr(jld2_dir::AbstractString; cfg=create_experiment_config(:full_model), reps::Int=3)
    _require_data(jld2_dir) || return nothing
    reps >= 1 || throw(ArgumentError("reps must be at least 1, got $reps"))

    println("Threads.nthreads() = ", Threads.nthreads())
    fields = load_greb_jld2!(jld2_dir; dataset=:ncep)

    # Warm-up, so compilation is not timed.
    redirect_stdout(devnull) do
        greb_model!(RunSpec(scnr=0), deepcopy(cfg); jld2_dir=jld2_dir, fields=deepcopy(fields))
    end

    times = Float64[]
    for r in 1:reps
        fields_r = deepcopy(fields)  # the model mutates fields
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

"Time each physics stage of one timestep; returns `(name, seconds per call)`."
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

"Time `year` in a subprocess per thread count; returns seconds per run for each."
function sweep_threads(jld2_dir::AbstractString; thread_counts=(1, 2, 3, 4), reps::Int=3)
    _require_data(jld2_dir) || return nothing

    script = @__FILE__

    results = Dict{Int,Vector{Float64}}()
    for n in thread_counts
        println("--- -t $n ---")
        cmd = `$JULIA_BIN --project=$REPO -t $n $script year $jld2_dir $reps`
        output = read(cmd, String)
        print(output)
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

"Bytes allocated by one `tendencies!` call, checked against the test budget."
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

const _MODES = ("year", "stages", "threads", "alloc")

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

    reps = length(rest) >= 2 ? parse_reps(rest[2]) : nothing

    if mode == "year"
        time_1yr(jld2_dir; reps=something(reps, 3))
    elseif mode == "stages"
        time_stages(jld2_dir; reps=something(reps, 2000))
    elseif mode == "threads"
        sweep_threads(jld2_dir; reps=something(reps, 3))
    elseif mode == "alloc"
        reps === nothing || error("the alloc mode takes no reps argument")
        check_allocations(jld2_dir)
    end
end
