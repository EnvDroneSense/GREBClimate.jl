# Profiles of GREBClimate.jl: where the time and the memory of a run go.
#
#   julia --project=. -t 1 benchmark/profile.jl [mode] [jld2_dir] [--years=N] [--samples=N] [--out=DIR] [--viewer=pprof]
#   julia --project=. benchmark/profile.jl compare <run before> <run after>
#
#   step     - sample control runs of `--years` simulated years (default 50), repeated
#              until `--samples` samples are collected (default 10000): the time loop
#   setup    - the same with 1-year runs: the cost of one `greb_model!` call
#   allocs   - record every allocation of one run of `--years` years (default 1)
#   dispatch - check one model step for runtime dispatch; static, nothing is timed
#   compare  - the change in shares between two step or setup runs; a run is a
#              report directory or its name under benchmark/profiles/
#
# Reports go to benchmark/profiles/<date>-<commit>-<mode>/ (or `--out`):
#
#   step, setup - header.txt, category.txt, owned.txt, flat.txt, tree.txt,
#                 stacks.folded, and profile.pb.gz with `--viewer=pprof`
#   allocs      - header.txt, allocs.txt
#   dispatch    - header.txt, dispatch.txt
#
# Run step and setup at `-t 1`: on Windows the sampler records the first thread
# only, so work on other threads shows up as waiting. `--viewer=pprof` needs
# PProf.jl and `dispatch` needs JET.jl in the default environment; neither is a
# dependency of the package.

using GREBClimate
using GREBClimate: tendencies!, time_loop!, init_model!
using Dates
using Printf
using Profile

include("common.jl")

const SRC_DIR = normpath(joinpath(REPO, "src"))
const PROFILES_DIR = joinpath(REPO, "benchmark", "profiles")
const WORKLOAD = "full_model control, stored corrections"

# Years per run and the most runs of each sampling mode. The sampler's rate
# varies, so a run is repeated until the sample target is reached.
const CPU_MODES = (step=(years=50, max_runs=8), setup=(years=1, max_runs=400))
const TARGET_SAMPLES = 10_000
const MIN_SAMPLES = 5_000        # below this only the largest rows mean anything
const PROFILE_BUFFER = 10_000_000
const META_WORDS = 6             # threadid, taskid, cycle clock, sleep state, 0, 0 close each sample

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
bump!(counts::AbstractDict, key, by::Int=1) = (counts[key] = get(counts, key, 0) + by; nothing)

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
        if data[i] == 0 && i - start + 1 >= META_WORDS && data[i-1] == 0
            push!(samples, (start:i-META_WORDS, i-META_WORDS+1:i))
            start = i + 1
        end
    end
    return samples
end

"""
Copy of `data` with the frames outside the package's outermost frame removed
from every sample, so a tree starts at the model and not at Julia's start-up.
A sample with no package frame is kept whole. `frames_at` maps an instruction
pointer to its frames.
"""
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

"""
The distinct stacks of a profile as `frames => samples` pairs: frames innermost
first, inlined frames included, cut at the package's outermost frame. A profile
of the time loop has far fewer distinct stacks than samples, so the reports are
built from these.
"""
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

"""
The innermost package frame of `stack` and the function it belongs to, or
`nothing` if the stack has no package frame.
"""
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

# ----- Reports of a sampled run -----

"`count` of `n` samples as a percentage and its standard error in percentage points, as padded strings."
function share(count::Int, n::Int)
    p = count / n
    return @sprintf("%7.2f", 100p), @sprintf("%6.2f", 100 * sqrt(p * (1 - p) / n))
end

"""
Write `owned.txt`: samples per package function and per package line, where a
sample belongs to its innermost package frame, so library time is charged to
the line that called it.
"""
function write_owned(path::AbstractString, stacks; mincount::Int, ms_per_run::Float64)
    n = sum(last, stacks)
    self = Dict{String,Int}()
    onstack = Dict{String,Int}()
    lines = Dict{String,Int}()
    for (stack, count) in stacks
        o = owner(stack)
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

    open(path, "w") do io
        println(io, "Samples: ", n, ". Shares in percent, +/- is one standard error in percentage points.")
        println(io, "Rows below ", mincount, " samples are left out of the tables but counted in the total.\n")
        println(io, "By function. self: the function holds the innermost package frame. on stack: it is anywhere in the stack.")
        println(io, "ms: the on-stack share of one run (", round(ms_per_run, digits=1), " ms).")
        println(io, "   self  share%    +/-  on stack  share%        ms  function")
        names = sort!(collect(union(keys(self), keys(onstack))); by=k -> (-get(self, k, 0), -get(onstack, k, 0), k))
        for name in names
            s, t = get(self, name, 0), get(onstack, name, 0)
            max(s, t) >= mincount || continue
            p, se = share(s, n)
            println(io, lpad(s, 7), " ", p, " ", se, lpad(t, 10), " ", share(t, n)[1],
                @sprintf("%10.2f", ms_per_run * t / n), "  ", name)
        end
        println(io, lpad(n, 7), " ", share(n, n)[1], "         total\n")

        println(io, "By line: the innermost package frame. A `@turbo` loop is one row, at its `@turbo for` line.")
        println(io, "   self  share%    +/-  line")
        for (line, count) in ranked(lines)
            count >= mincount || break
            p, se = share(count, n)
            println(io, lpad(count, 7), " ", p, " ", se, "  ", line)
        end
    end
end

"""
Write `category.txt`: samples per kind of cost, judged by the innermost frame,
with the largest frames and the package functions they are charged to.
"""
function write_category(path::AbstractString, stacks)
    n = sum(last, stacks)
    classes = Dict{String,Int}()
    frames = Dict{String,Dict{String,Int}}()
    owners = Dict{String,Dict{String,Int}}()
    for (stack, count) in stacks
        isempty(stack) && continue
        fr = first(stack)
        class = classify(fr)
        o = owner(stack)
        bump!(classes, class, count)
        bump!(get!(Dict{String,Int}, frames, class), label(fr), count)
        bump!(get!(Dict{String,Int}, owners, class), o === nothing ? NO_OWNER : o[2], count)
    end

    open(path, "w") do io
        println(io, "Samples: ", n, ". Each sample is classed by its innermost frame.")
        println(io, "Shares in percent, +/- is one standard error in percentage points.\n")
        println(io, "samples  share%    +/-  class")
        for (class, count) in ranked(classes)
            p, se = share(count, n)
            println(io, lpad(count, 7), " ", p, " ", se, "  ", class)
        end
        total = sum(values(classes))
        println(io, lpad(total, 7), " ", share(total, n)[1], "         total")
        for (class, _) in ranked(classes)
            println(io, "\n", class)
            println(io, "  frames:     ", largest(frames[class]))
            println(io, "  charged to: ", largest(owners[class]))
        end
    end
end

"""
Write `stacks.folded`: one line per distinct stack, outermost frame first, the
frames separated by `;` and followed by the sample count. This is the collapsed
format flame graph tools read (speedscope, FlameGraph, inferno).
"""
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

"""
An optional package, loaded from the active load path. PProf.jl and JET.jl are
not dependencies of this script, so they have to be installed in the default
environment.
"""
function load_optional(name::AbstractString, uuid::AbstractString)
    try
        Base.require(Base.PkgId(Base.UUID(uuid), name))
    catch
        error("this needs $name.jl in the default environment: julia -e 'using Pkg; Pkg.add(\"$name\")'")
    end
end

load_pprof() = load_optional("PProf", "e4faabce-9ead-11e9-39d9-4379958e3056")
load_jet() = load_optional("JET", "c3a54625-cd67-489e-a8e7-0a5a0ff4e31b")

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

"One aligned `name: value` line of a header."
header_line(name::AbstractString, value) = string(rpad(name * ":", 13), value, "\n")

"""
Make the directory for a run's reports, `out_dir` or a new one under
`benchmark/profiles/` (numbered if the name is taken), and return it with the
lines every header starts with.
"""
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
    header = header_line("mode", "$mode ($description)") *
        header_line("commit", commit * (dirty ? " (tracked files modified)" : "")) *
        header_line("date", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM")) *
        header_line("julia", VERSION)
    return dir, header
end

"Write `header` to the run's `header.txt`, print it, and return `dir`."
function close_report(dir::AbstractString, header::AbstractString)
    write(joinpath(dir, "header.txt"), header)
    print(header)
    println("reports written to ", dir)
    return dir
end

"Run `f` with the run's progress lines and log messages discarded; they are not part of the workload."
quiet(f) = redirect_stdout(devnull) do
    Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
end

"""
The configuration and the loaded fields that every mode runs: `:full_model` on
the stored flux corrections, no spin-up. One year is run first, so compilation
is not profiled. `greb_model!` restores what it overwrites in `fields`, so the
same fields serve every run; a copy per run would leave garbage for the
collector to clear inside the profiled run.
"""
function workload(jld2_dir::AbstractString)
    isdir(jld2_dir) || error("JLD2 data directory not found: $jld2_dir. Set GREB_DATA or pass a path.")
    cfg = preset(:full_model; corrections=Stored())
    fields = load_climatology(jld2_dir; dataset=:ncep)
    quiet() do
        greb_model!(RunSpec(ctrl=1, scnr=0), cfg; jld2_dir=jld2_dir, fields=fields)
    end
    GC.gc()
    return cfg, fields
end

"""
Sample `years`-year control runs, repeated until `target` samples are collected
or the mode's run limit is reached, and write the reports to `out_dir`; returns
it. `mode` is `:step` (long runs, so the time loop dominates) or `:setup`
(1-year runs, so the cost of one `greb_model!` call shows).
"""
function profile_cpu(jld2_dir::AbstractString, mode::Symbol; years::Int=CPU_MODES[mode].years,
        target::Int=TARGET_SAMPLES, out_dir=nothing, viewer=nothing)
    years >= 1 || throw(ArgumentError("years must be at least 1, got $years"))
    viewer in (nothing, "pprof") || throw(ArgumentError("unknown viewer $(repr(viewer)); expected pprof"))
    PProf = viewer === nothing ? nothing : load_pprof()   # before the run, so a missing package fails early
    if Threads.nthreads() > 1
        @warn "Profiling with $(Threads.nthreads()) threads: on Windows only the first thread is sampled. Use -t 1."
    end
    cfg, fields = workload(jld2_dir)

    Profile.clear()
    Profile.init(n=PROFILE_BUFFER, delay=0.001)
    elapsed = 0.0
    runs = 0
    while runs < CPU_MODES[mode].max_runs &&
            (runs == 0 || length(split_samples(Profile.fetch(include_meta=true))) < target)
        elapsed += @elapsed quiet() do
            @profile greb_model!(RunSpec(ctrl=years, scnr=0), cfg; jld2_dir=jld2_dir, fields=fields)
        end
        runs += 1
        Profile.is_buffer_full() && error("the profile buffer filled up; lower --years or --samples")
    end

    data, lidict = Profile.retrieve(include_meta=true)
    samples = split_samples(data)
    n = length(samples)
    n > 0 || error("the profiler recorded no samples")
    threads_seen = sort!(unique(Int(data[first(meta)]) for (_, meta) in samples))
    ms_per_run = 1000 * elapsed / runs

    dir, header = open_report(string(mode), WORKLOAD, out_dir)
    header *= header_line("threads", "$(Threads.nthreads()) (sampled: $(join(threads_seen, ", ")))") *
        header_line("runs", "$runs of $years simulated year(s)") *
        header_line("elapsed", "$(round(elapsed, digits=2)) s ($(round(ms_per_run, digits=1)) ms per run)") *
        header_line("samples", n) *
        header_line("interval", "$(round(1000 * elapsed / n, digits=2)) ms per sample")

    frames_at = Dict(ip => Frame.(frames) for (ip, frames) in lidict)
    trimmed = trim_to_package(data, samples, frames_at)
    stacks = distinct_stacks(data, samples, frames_at)
    mincount = max(2, n ÷ 1000)   # rows below 0.1 percent are noise
    write_category(joinpath(dir, "category.txt"), stacks)
    write_owned(joinpath(dir, "owned.txt"), stacks; mincount, ms_per_run)
    # Sorted by self time, largest last.
    write_report(joinpath(dir, "flat.txt"), trimmed, lidict; format=:flat, sortedby=:overhead, mincount)
    write_report(joinpath(dir, "tree.txt"), trimmed, lidict; format=:tree, mincount, noisefloor=2)
    write_folded(joinpath(dir, "stacks.folded"), stacks)
    PProf === nothing || Base.invokelatest(PProf.pprof, trimmed, lidict; web=false,
        out=joinpath(dir, "profile.pb.gz"), full_signatures=false)

    close_report(dir, header)
    n >= MIN_SAMPLES ||
        @warn "Only $n samples (fewer than $MIN_SAMPLES): shares of small rows are not reliable. Raise --samples or --years."
    return dir
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

"`bytes` with a unit, for a table."
function format_bytes(bytes::Integer)
    bytes >= 2^20 ? @sprintf("%9.2f MB", bytes / 2^20) : @sprintf("%9.2f kB", bytes / 2^10)
end

"""
Record every allocation of one `years`-year control run and write `allocs.txt`
to `out_dir`: bytes and counts per package function and per package line, an
allocation belonging to the innermost package frame of its stack. Returns
`out_dir`.
"""
function profile_allocs(jld2_dir::AbstractString; years::Int=1, out_dir=nothing)
    years >= 1 || throw(ArgumentError("years must be at least 1, got $years"))
    cfg, fields = workload(jld2_dir)

    Profile.Allocs.clear()
    quiet() do
        Profile.Allocs.@profile sample_rate = 1 greb_model!(RunSpec(ctrl=years, scnr=0), cfg;
            jld2_dir=jld2_dir, fields=fields)
    end
    allocs = Profile.Allocs.fetch().allocs
    isempty(allocs) && error("the allocation profiler recorded nothing")

    frame_of = Dict{Base.StackTraces.StackFrame,Frame}()
    by_func = Dict{String,AllocRow}()
    by_line = Dict{String,AllocRow}()
    total, in_loop = AllocRow(), AllocRow()
    for a in allocs
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
    amount(row) = string(row.count, ", ", strip(format_bytes(row.bytes)))
    largest_first(rows) = sort!(collect(rows); by=kv -> (-last(kv).bytes, first(kv)))

    dir, header = open_report("allocs", WORKLOAD, out_dir)
    header *= header_line("runs", "1 of $years simulated year(s)") *
        header_line("allocations", amount(total)) *
        header_line("in loop", amount(in_loop) * " (with `time_loop!` on the stack)")

    open(joinpath(dir, "allocs.txt"), "w") do io
        println(io, "Every allocation of one run, charged to the innermost package frame of its stack.")
        println(io, "Total: ", amount(total), ". With `time_loop!` on the stack: ", amount(in_loop), ".\n")
        println(io, "By function")
        println(io, "       bytes  share%   count  function")
        for (name, row) in largest_first(by_func)
            println(io, format_bytes(row.bytes), " ", @sprintf("%7.2f", 100 * row.bytes / total.bytes),
                lpad(row.count, 8), "  ", name)
        end
        println(io, "\nBy line (the 40 largest)")
        println(io, "       bytes  share%   count  line, and the type holding the most bytes there (bytes)")
        for (line, row) in first(largest_first(by_line), 40)
            println(io, format_bytes(row.bytes), " ", @sprintf("%7.2f", 100 * row.bytes / total.bytes),
                lpad(row.count, 8), "  ", line, "  [", first(largest(row.types, 1), 80), "]")
        end
    end
    return close_report(dir, header)
end

"""
Check one model step for runtime dispatch with JET's optimization analysis and
write `dispatch.txt` to `out_dir`; returns it. This is a static check of the
compiled code: nothing is timed, and the serial and the threaded branch of
`tendencies!` are both covered. A report names a call whose target is not
known at compile time, inside the package only.
"""
function profile_dispatch(jld2_dir::AbstractString; out_dir=nothing)
    JET = load_jet()   # before the set-up, so a missing package fails early
    cfg, fields = workload(jld2_dir)
    r = resolve(cfg; jld2_dir)
    CO2 = init_model!(r, fields).CO2_ctrl
    state = ModelState()
    ws = ModelWorkspace()
    acc = GREBClimate.MonthlyAccumulator()
    timestate = TimeState(1, 1)
    Ts = fields.Ts_clim[:, :, 1]
    Ta = copy(Ts)
    To = fields.To_clim[:, :, 1]
    q = fields.q_clim[:, :, 1]
    records = GREBClimate.MonthlyRecord[]

    step = () -> GREBClimate.time_loop!(1, 1, CO2, 1, 0, Ts, Ta, q, To, records, fields, state, ws, acc, timestate, r)
    result = Base.invokelatest(JET.report_opt, step, (); target_modules=(GREBClimate,))
    n = length(Base.invokelatest(JET.get_reports, result))

    dir, header = open_report("dispatch", "one time_loop! step of full_model", out_dir)
    open(joinpath(dir, "dispatch.txt"), "w") do io
        println(io, "JET optimization analysis of one `time_loop!` step, reports from the package only.")
        println(io, "A report is a call that is dispatched at run time. Reports: ", n, ".\n")
        n > 0 && Base.invokelatest(show, IOContext(io, :color => false), MIME"text/plain"(), result)
    end
    return close_report(dir, header * header_line("reports", n))
end

# ----- Comparing two runs -----

"A run's report directory: `path` itself, or the folder of that name under `benchmark/profiles/`."
function run_dir(path::AbstractString)
    isdir(path) && return path
    inside = joinpath(PROFILES_DIR, path)
    isdir(inside) || error("no report directory at $path or $inside")
    return inside
end

"""
The tables of a `step` or `setup` run read back from its reports:
`(commit, samples, functions, lines, classes)`, the last three mapping a row
name to its self samples.
"""
function read_run(dir::AbstractString)
    header = read(joinpath(dir, "header.txt"), String)
    commit = strip(match(r"commit:\s*(.*)", header).captures[1])
    samples = parse(Int, match(r"samples:\s*(\d+)", header).captures[1])

    functions, lines, classes = Dict{String,Int}(), Dict{String,Int}(), Dict{String,Int}()
    table = functions
    for row in eachline(joinpath(dir, "owned.txt"))
        startswith(row, "By line") && (table = lines)
        m = table === functions ?
            match(r"^\s*(\d+)\s+[\d.]+\s+[\d.]+\s+\d+\s+[\d.]+(?:\s+[\d.]+)?\s+([^\d\s].*)$", row) :
            match(r"^\s*(\d+)\s+[\d.]+\s+[\d.]+\s+(\S+:\d+\s+\S.*)$", row)
        m === nothing || (table[m.captures[2]] = parse(Int, m.captures[1]))
    end
    for row in eachline(joinpath(dir, "category.txt"))
        isempty(strip(row)) && !isempty(classes) && break   # the class table ends at the first blank line
        m = match(r"^\s*(\d+)\s+[\d.]+\s+[\d.]+\s+([^\d\s].*)$", row)
        m === nothing || (classes[m.captures[2]] = parse(Int, m.captures[1]))
    end
    return (; commit, samples, functions, lines, classes)
end

"""
Print one table of `compare`: the share of each row in runs `a` and `b`, the
difference in percentage points and its standard error. A `*` marks a
difference of three standard errors or more, a `?` a row one run does not list.
"""
function print_changes(io::IO, title::AbstractString, a::Dict{String,Int}, na::Int, b::Dict{String,Int}, nb::Int;
        top::Int=25)
    rows = map(collect(union(keys(a), keys(b)))) do name
        pa, pb = get(a, name, 0) / na, get(b, name, 0) / nb
        se = sqrt(pa * (1 - pa) / na + pb * (1 - pb) / nb)
        (; name, pa, pb, diff=pb - pa, se)
    end
    sort!(rows; by=r -> (-abs(r.diff), r.name))
    println(io, title)
    println(io, "      A%      B%    diff    +/-     row")
    for r in first(rows, top)
        listed = haskey(a, r.name) && haskey(b, r.name)
        flag = !listed ? " ? " : abs(r.diff) >= 3 * r.se ? " * " : "   "
        println(io, @sprintf("%8.2f%8.2f%+8.2f%7.2f", 100r.pa, 100r.pb, 100r.diff, 100r.se), flag, " ", r.name)
    end
    println(io)
end

"""
Compare the shares of two `step` or `setup` runs, `a` before and `b` after,
largest change first. Shares, not times: a row's share also falls when
another row grows.
"""
function compare_runs(dir_a::AbstractString, dir_b::AbstractString; io::IO=stdout)
    a, b = read_run(run_dir(dir_a)), read_run(run_dir(dir_b))
    println(io, "A: ", dir_a, "  commit ", a.commit, ", ", a.samples, " samples")
    println(io, "B: ", dir_b, "  commit ", b.commit, ", ", b.samples, " samples")
    println(io, "Self shares in percent; diff is B minus A in percentage points, +/- its standard error.")
    println(io, "* marks a difference of three standard errors or more. ? marks a row that one run's report")
    println(io, "does not list (below its cut-off, or a line that moved); it counts as zero there.\n")
    print_changes(io, "By class", a.classes, a.samples, b.classes, b.samples)
    print_changes(io, "By function", a.functions, a.samples, b.functions, b.samples)
    print_changes(io, "By line", a.lines, a.samples, b.lines, b.samples)
    return nothing
end

# ----- Command line -----

const _PROFILE_MODES = ("step", "setup", "allocs", "dispatch", "compare")

# The flags each mode takes.
const _MODE_FLAGS = (step=("years", "samples", "out", "viewer"), setup=("years", "samples", "out", "viewer"),
    allocs=("years", "out"), dispatch=("out",), compare=())

function main(args::Vector{String})
    mode, rest = if isempty(args) || startswith(args[1], "--")
        ("step", args)
    elseif args[1] in _PROFILE_MODES
        (args[1], args[2:end])
    else
        error("unknown mode $(repr(args[1])); expected one of $(join(_PROFILE_MODES, ", "))")
    end

    flags, rest = split_flags(rest)
    for name in keys(flags)
        name in _MODE_FLAGS[Symbol(mode)] || error("the $mode mode takes no --$name flag")
    end
    if mode == "compare"
        length(rest) == 2 || error("compare takes two report directories: the run before and the run after")
        return compare_runs(rest[1], rest[2])
    end

    jld2_dir = !isempty(rest) ? rest[1] : default_data_dir()
    out_dir = get(flags, "out", nothing)
    years = haskey(flags, "years") ? parse_nonneg_int("years", flags["years"]) : nothing
    if mode == "allocs"
        profile_allocs(jld2_dir; years=something(years, 1), out_dir)
    elseif mode == "dispatch"
        profile_dispatch(jld2_dir; out_dir)
    else
        target = haskey(flags, "samples") ? parse_nonneg_int("samples", flags["samples"]) : TARGET_SAMPLES
        profile_cpu(jld2_dir, Symbol(mode); years=something(years, CPU_MODES[Symbol(mode)].years), target, out_dir,
            viewer=get(flags, "viewer", nothing))
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
