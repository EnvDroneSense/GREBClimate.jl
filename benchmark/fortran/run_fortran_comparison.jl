# Compare GREBClimate.jl with the original Fortran GREB on the same workload.
#
#   julia --project=benchmark benchmark/fortran/run_fortran_comparison.jl [mode] [--reps=N]
#
#   compare  - time each Fortran build and Julia, and report the Julia peak memory
#              (default, --reps=3)
#   verify   - check both models give the same result
#   io       - time the Fortran's output writes
#   build    - compile the Fortran builds
#   driver   - the Julia run, started by the other modes in a child process
#
# Needs gfortran and the Fortran source tree (GREB_GFORTRAN, GREB_FORTRAN_DIR).
# Exits 1 when a mode fails: missing tooling, a failed gate or self-check.
# Peak memory of Julia is its own Sys.maxrss (the working set, with runtime, GC heap and packages).

using GREBClimate: RunSpec, SpinUp, greb_data_dir, greb_model!, load_climatology, preset
using Random

# The child process of the timed Julia run ("driver") loads the package only: Bench would add its load time.
const DRIVER = !isempty(ARGS) && ARGS[1] == "driver"
if !DRIVER
    include(joinpath(@__DIR__, "..", "..", "tools", "ToolBox.jl"))
    using .ToolBox
    include(joinpath(@__DIR__, "..", "Bench.jl"))
    using Statistics: median
end
const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const JULIA_BIN = joinpath(Sys.BINDIR, Base.julia_exename())

# Fortran grid and output layout.
const XDIM = 96
const YDIM = 48
const FIELDS_PER_RECORD = 8   # Tmm, Tamm, Tomm, qmm, icmm, prmm, evamm, qcrclmm

# Workload: the default of `run.greb.scenarios.csh`.
const TIME_FLUX = 3
const TIME_CTRL = 3
const TIME_SCNR = 50
# Fortran EXP and its Julia counterpart; GREB_FORTRAN_EXP picks another pair.
const EXPERIMENTS = Dict(20 => :co2_double, 40 => :regional_co2_nh, 41 => :regional_co2_sh,
    42 => :regional_co2_tropics, 43 => :regional_co2_extratropics)
const LOG_EXP = parse(Int, get(ENV, "GREB_FORTRAN_EXP", "20"))
const JULIA_EXPERIMENT = EXPERIMENTS[LOG_EXP]
const THREADS = 2                     # OMP_NUM_THREADS / julia -t

const DATASET = :era        # the Fortran's climatology (log_clim = 0)
const EVAPORATION = :skin   # the Fortran's evaporation scheme (log_eva = 0)

# `verify` tolerances on the scenario Ts anomaly.
const GATE_RELATIVE_TOLERANCE = 0.01   # final-year global mean, relative
const GATE_MEAN_ABS_K = 0.005          # pointwise mean |Fortran - Julia|, all months
const GATE_MAX_ABS_K = 0.2             # pointwise max  |Fortran - Julia|, all months
const GATE_MIN_SIGNAL_K = 0.1          # the Fortran anomaly must be at least this, or two near-zero runs would pass

# The Windows stack is too small otherwise; this is a PE linker option.
const STACK_FLAG = Sys.iswindows() ? ["-Wl,--stack,536870912"] : String[]

# `control` is a copy of `official`, used as a noise check.
const BUILDS = (
    official=["-Ofast", "-ffast-math", "-funroll-loops", "-fopenmp"],
    native=["-Ofast", "-ffast-math", "-funroll-loops", "-fopenmp", "-march=native"],
    control=["-Ofast", "-ffast-math", "-funroll-loops", "-fopenmp"]
)

# ── Toolchain and sources ───────────────────────────────────────────────────

"Path to gfortran, or `nothing`."
function find_gfortran()
    candidates = filter(!isempty, [get(ENV, "GREB_GFORTRAN", ""), "gfortran",
        joinpath(homedir(), "scoop", "apps", "mingw", "current", "bin", "gfortran.exe")])
    for c in candidates
        path = Sys.which(c)   # accepts a name or a path
        path !== nothing && success(ignorestatus(`$path --version`)) && return abspath(path)
    end
    nothing
end

"Path to the Fortran source tree, or `nothing`."
function fortran_source_dir()
    candidates = filter(!isempty, [get(ENV, "GREB_FORTRAN_DIR", ""),
        (normpath(joinpath(REPO, up, "greb-official-official")) for up in ("..", "../..", "../data"))...])
    for c in candidates
        isfile(joinpath(c, "greb.model.mscm.f90")) && isfile(joinpath(c, "greb.shell.mscm.f90")) &&
            isdir(joinpath(c, "input")) && return c
    end
    nothing
end

"Compiler, Fortran tree and dataset paths; an error names what is missing."
function require_tooling()
    t = (gfortran=find_gfortran(), fortran_dir=fortran_source_dir(), jld2_dir=data_dir())
    missing_parts = String[]
    t.gfortran === nothing && push!(missing_parts, "gfortran not found. Set GREB_GFORTRAN, or `scoop install mingw`.")
    t.fortran_dir === nothing && push!(missing_parts, "Official Fortran GREB tree not found. Set GREB_FORTRAN_DIR " *
                                                      "to a directory holding greb.model.mscm.f90, greb.shell.mscm.f90 and input/.")
    isdir(t.jld2_dir) || push!(missing_parts, "JLD2 data directory not found: $(t.jld2_dir). Set GREB_DATA.")
    isempty(missing_parts) || error(join(missing_parts, "\n"))
    t
end

# ── Build tree ──────────────────────────────────────────────────────────────

"Scratch directory for builds and runs."
bench_root() = normpath(abspath(get(ENV, "GREB_BENCH_DIR", joinpath(tempdir(), "greb_fortran_bench"))))

"Copy the Fortran sources to the scratch directory and link its input data."
function prepare_tree(fortran_dir::AbstractString)
    root = mkpath(bench_root())
    for f in ("greb.model.mscm.f90", "greb.shell.mscm.f90")
        cp(joinpath(fortran_dir, f), joinpath(root, f); force=true)
    end
    link = joinpath(root, "input")
    target = normpath(joinpath(fortran_dir, "input"))
    # Recreate the link when the source tree has changed, so the old input data is not used.
    if !isdir(link) || !endswith(normpath(readlink(link)), target)
        isdir(link) && rm(link)   # not recursive: removes the link, not its target
        Sys.iswindows() ? run(`cmd /c mklink /J $link $target`) : symlink(target, link; dir_target=true)
    end
    root
end

"Run `cmd` with its output in `log`; on failure print the log's tail, then rethrow."
function run_logged(cmd::Base.AbstractCmd, log::AbstractString)
    try
        run(pipeline(cmd; stdout=log, stderr=log))
    catch
        lines = isfile(log) ? readlines(log) : String[]
        println("\n  command failed; last lines of ", log, ":")
        foreach(l -> println("    ", l), lines[max(1, end-19):end])
        rethrow()
    end
end

"Compile every build; returns name => executable."
function build_all(gfortran::AbstractString, root::AbstractString)
    exes = Dict{Symbol,String}()
    for (name, flags) in pairs(BUILDS)
        moddir = mkpath(joinpath(root, "mod_$name"))
        exe = joinpath(root, "greb_$name.exe")
        srcs = joinpath.(root, ["greb.model.mscm.f90", "greb.shell.mscm.f90"])
        print("  building $name ... ")
        t = @elapsed run_logged(`$gfortran $flags $STACK_FLAG -J$moddir $srcs -o $exe`,
            joinpath(root, "build_$name.log"))
        println(round(t, digits=1), "s")
        exes[name] = exe
    end
    exes
end

# ── Running each side ───────────────────────────────────────────────────────

"Write the Fortran namelist and placeholder input files."
function write_namelist(workdir::AbstractString)
    write(joinpath(workdir, "namelist"), """
        &NUMERICS
        time_flux = $TIME_FLUX
        time_ctrl = $TIME_CTRL
        time_scnr = $TIME_SCNR
        /
        &PHYSICS
         log_exp = $LOG_EXP
         dradius = 0
        /
        """)
    foreach(f -> touch(joinpath(workdir, f)), ("nosolfile", "noco2file"))
end

"Run one Fortran build; returns wall-clock seconds."
function run_fortran(exe::AbstractString, workdir::AbstractString, compiler_bindir::AbstractString)
    mkpath(workdir)
    write_namelist(workdir)
    cmd = addenv(Cmd(`$exe`; dir=workdir), "OMP_NUM_THREADS" => string(THREADS),
        "PATH" => string(compiler_bindir, Sys.iswindows() ? ";" : ":", get(ENV, "PATH", "")))
    @elapsed run_logged(cmd, joinpath(workdir, "run.log"))
end

"The Julia side of the workload; with `dump_path` it also saves the scenario Ts."
function julia_run(jld2_dir::AbstractString, dump_path::AbstractString="")
    fields = load_climatology(jld2_dir; dataset=DATASET)
    cfg = preset(JULIA_EXPERIMENT; hydrology=(evaporation=EVAPORATION,), corrections=SpinUp(TIME_FLUX))
    result = greb_model!(RunSpec(ctrl=TIME_CTRL, scnr=TIME_SCNR), cfg; jld2_dir, fields)
    isempty(dump_path) || open(dump_path, "w") do out
        foreach(rec -> write(out, rec.Ts), result.scnr)
    end
    println("peak_rss_bytes ", Sys.maxrss())   # read by `julia_peak_mib`
end

"The peak memory in MiB that the driver wrote to `log`, or `nothing`."
function julia_peak_mib(log::AbstractString)
    m = match(r"^peak_rss_bytes (\d+)$"m, read(log, String))
    m === nothing ? nothing : round(parse(Int, m[1]) / 1024^2, digits=1)
end

"Arguments of the Julia process that runs `julia_run` (this script in `driver` mode)."
julia_args(jld2_dir, dump_path="") = ["--project=$REPO", "-t", string(THREADS), @__FILE__, "driver", jld2_dir, (isempty(dump_path) ? [] : [dump_path])...]

"Run `julia_run` in a new process, with its output in `log`; returns wall-clock seconds."
run_julia(jld2_dir::AbstractString, log::AbstractString; dump_path::AbstractString="") = @elapsed run_logged(`$JULIA_BIN $(julia_args(jld2_dir, dump_path))`, log)

# ── Correctness gate ────────────────────────────────────────────────────────

"Read the scenario Ts fields from the Fortran's `scenario.bin`."
function read_fortran_scenario_ts(path::AbstractString, nmonths::Int)
    reclen = XDIM * YDIM * sizeof(Float32)
    out = Array{Float32,3}(undef, XDIM, YDIM, nmonths)
    open(path) do io
        for i in 1:nmonths
            seek(io, (FIELDS_PER_RECORD * i - 8) * reclen)   # Ts is the first field of month i
            read!(io, view(out, :, :, i))
        end
    end
    out
end

"`cos(latitude)` weight per latitude row."
area_weights() = [cosd(-90 + (j - 0.5) * 180 / YDIM) for j in 1:YDIM]

"Area-weighted global mean of the last 12 months."
function final_year_gmean(field::Array{Float32,3})
    w = area_weights()
    n = size(field, 3)
    sum(Float64(field[i, j, m]) * w[j] for m in (n-11):n, j in 1:YDIM, i in 1:XDIM) / (12 * XDIM * sum(w))
end

"Run both models once and check they agree; `exes` reuses built Fortran executables. Returns the results and `pass`."
function verify(; exes=nothing)
    t = require_tooling()
    root = prepare_tree(t.fortran_dir)
    exes === nothing && (exes = build_all(t.gfortran, root))
    nmonths = 12 * TIME_SCNR

    println("  running Fortran (untimed) ...")
    workdir = joinpath(root, "verify_work")
    run_fortran(exes[:official], workdir, dirname(t.gfortran))
    println("  running Julia (untimed) ...")
    dump = joinpath(root, "julia_scnr_ts.bin")
    rm(dump; force=true)
    run_julia(t.jld2_dir, joinpath(root, "julia_verify.log"); dump_path=dump)

    ftn = read_fortran_scenario_ts(joinpath(workdir, "scenario.bin"), nmonths)
    jl = read!(dump, Array{Float32,3}(undef, XDIM, YDIM, nmonths))
    fg, jg = final_year_gmean(ftn), final_year_gmean(jl)
    diff = abs.(ftn .- jl)
    rel = abs(fg - jg) / max(abs(fg), eps())
    meandiff, maxdiff = sum(diff) / length(diff), maximum(diff)
    # (what is checked, value, limit, passed)
    checks = [("Fortran final-year anomaly (K, at least)", abs(fg), GATE_MIN_SIGNAL_K, abs(fg) > GATE_MIN_SIGNAL_K),
        ("final-year mean, relative difference", rel, GATE_RELATIVE_TOLERANCE, rel <= GATE_RELATIVE_TOLERANCE),
        ("mean |Fortran - Julia| over all months (K)", meandiff, GATE_MEAN_ABS_K, meandiff <= GATE_MEAN_ABS_K),
        ("max  |Fortran - Julia| over all months (K)", maxdiff, GATE_MAX_ABS_K, maxdiff <= GATE_MAX_ABS_K)]
    println("
Correctness gate - did both models run the same experiment?")
    println("  final-year global-mean Ts anomaly: Fortran ", round(fg, digits=4), " K, Julia ", round(jg, digits=4), " K")
    for (what, value, limit, ok) in checks
        println("  ", rpad(what, 45), rpad(round(value, sigdigits=3), 9), "limit ", rpad(limit, 6), ok ? "PASS" : "FAIL")
    end
    pass = all(last, checks)
    println("  overall -> ", pass ? "PASS" : "FAIL")
    (; pass, fortran=fg, julia=jg)
end

# ── I/O asymmetry ───────────────────────────────────────────────────────────

"Time writing the Fortran's output volume; returns seconds."
function io_cost()
    nrec = FIELDS_PER_RECORD * 12 * TIME_SCNR
    buf = rand(Float32, XDIM, YDIM)
    path = joinpath(mkpath(bench_root()), "io_probe.bin")
    t = @elapsed open(path, "w") do io
        foreach(_ -> write(io, buf), 1:nrec)
        flush(io)
    end
    rm(path; force=true)
    mb = nrec * sizeof(buf) / 1024^2
    println("I/O the Fortran does and Julia does not:")
    println("  ", nrec, " records, ", round(mb, digits=1), " MB, ", round(t, digits=2), "s  (",
        round(mb / t, digits=0), " MB/s)")
    t
end

# ── The comparison ──────────────────────────────────────────────────────────

"Time each Fortran build and Julia over `reps` shuffled trials; returns whether the gate and the self-check passed."
function compare(; reps::Int=3)
    t = require_tooling()
    root = prepare_tree(t.fortran_dir)
    println("Build tree: ", root)
    exes = build_all(t.gfortran, root)
    verify(; exes).pass || (println("
Correctness gate failed - not reporting timings."); return false)

    println("
Workload: EXP=", LOG_EXP, " / ", JULIA_EXPERIMENT, ", flux=", TIME_FLUX, " ctrl=", TIME_CTRL,
        " scnr=", TIME_SCNR, " (", TIME_FLUX + TIME_CTRL + TIME_SCNR, " simulated years)")
    println("Threads: OMP_NUM_THREADS=", THREADS, " (Fortran) / -t ", THREADS, " (Julia)")
    for (name, flags) in pairs(BUILDS)
        println("Flags [", name, "]: ", join([flags; STACK_FLAG], " "))
    end
    println()

    variants = [collect(keys(BUILDS)); :julia]
    times = Dict(v => Float64[] for v in variants)
    peaks = Float64[]   # MiB of Julia
    before = Bench.machine_header()
    for r in 1:reps
        println("--- trial $r/$reps ---")
        for v in shuffle(variants)
            s = v === :julia ? run_julia(t.jld2_dir, joinpath(root, "julia_timed.log")) :
                run_fortran(exes[v], joinpath(root, "work_$v"), dirname(t.gfortran))
            push!(times[v], s)
            println("  ", rpad(v, 9), round(s, digits=2), "s")
            v === :julia && push!(peaks, something(julia_peak_mib(joinpath(root, "julia_timed.log")), NaN))
        end
    end

    r2(x) = round(x, digits=2)
    println("
Results (", reps, " reps):")
    for v in variants
        x = times[v]
        println("  ", rpad(v, 9), "min ", rpad(r2(minimum(x)), 7), "median ", rpad(r2(median(x)), 7), "max ", r2(maximum(x)))
    end
    best = Dict(v => minimum(times[v]) for v in variants)
    println("
Fortran : Julia (min of trials)")
    println("  official / julia = ", r2(best[:official] / best[:julia]), "x")
    println("  native   / julia = ", r2(best[:native] / best[:julia]), "x")
    selfcheck = best[:control] / best[:official]
    selfcheck_ok = abs(selfcheck - 1) <= 0.05
    println("
Self-check (control / official, identical builds): ", round(selfcheck, digits=3), "x  ",
        selfcheck_ok ? "ok" : "WARNING: not within 5% of 1.00x - DISCARD THIS RUN.")
    println("
Julia peak memory (Sys.maxrss, max of trials): ", maximum(peaks), " MiB, runtime and packages included.")
    println("The Julia time is a fresh process: start-up, package load and compilation are inside it.")
    println("The Fortran also writes ~88 MB that Julia does not; `io` mode bounds its cost.
")
    quiet = Bench.machine_footer(before, times[:official])
    for v in variants
        Bench.record_result("fortran", string(v), times[v], quiet; selfcheck=v === :control ? selfcheck : nothing, threads=THREADS)
    end
    selfcheck_ok
end

const FLAGS = Dict("compare" => ("reps",), "verify" => (), "io" => (), "build" => ())

"Run the mode named by `args`; returns the exit code."
function main(args::Vector{String})
    if !isempty(args) && args[1] == "driver"
        julia_run(args[2:end]...)
        return 0
    end
    mode, o, positional = parse_args(args, FLAGS, "compare")
    isempty(positional) || error("the $mode mode takes no positional argument, got: ", join(positional, " "))
    ok = if mode == "compare"
        compare(; reps=intflag(o, "reps", 3; min=1))
    elseif mode == "verify"
        verify().pass
    elseif mode == "io"
        io_cost() >= 0
    else  # build
        t = require_tooling()
        build_all(t.gfortran, prepare_tree(t.fortran_dir))
        true
    end
    return ok ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
