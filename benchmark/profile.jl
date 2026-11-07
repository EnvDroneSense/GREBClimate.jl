# Profiles of GREBClimate.jl: where the time and the memory of a run go.
#
#   julia --project=. -t 1 benchmark/profile.jl [mode] [jld2_dir] [--years=N] [--samples=N] [--out=DIR] [--viewer=pprof]
#   julia --project=. benchmark/profile.jl compare <run before> <run after>
#
#   step     - sample control runs of `--years` simulated years (default 25), repeated
#              until `--samples` samples are collected (default 10000): the time loop
#   setup    - the same with 1-year runs: the cost of one `greb_model!` call
#   allocs   - record every allocation of one run of `--years` years (default 1)
#   dispatch - check one model step for runtime dispatch; static, nothing is timed
#   compare  - the change in shares between two step or setup runs; a run is a
#              report directory or its name under benchmark/profiles/
#
# Reports go to benchmark/profiles/<date>-<commit>-<mode>/ (or `--out`): report.md
# for every mode; step and setup add counts.tsv, flat.txt, tree.txt, stacks.folded
# and, with `--viewer=pprof`, profile.pb.gz. Details: benchmark/README.md.
#
# Run step and setup at `-t 1` (the Windows sampler records the first thread
# only). `--viewer=pprof` needs PProf.jl and `dispatch` needs JET.jl in the
# default environment. Load one harness per Julia session: this file and
# run_benchmarks.jl both define `main`.

using GREBClimate
using Dates
using Printf
using Profile

include("common.jl")

const SRC_DIR = normpath(joinpath(REPO, "src"))
const PROFILES_DIR = joinpath(REPO, "benchmark", "profiles")
const WORKLOAD = "full_model control, stored corrections"

# Years per run and the most runs of each sampling mode. The sampler's rate
# varies, so a run is repeated until the sample target is reached.
const CPU_MODES = (step=(years=25, max_runs=8), setup=(years=1, max_runs=400))
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

"`lidict` with the inlined frames outside the outermost package frame cut, as `trim_to_package` does for the samples."
function trim_lidict(lidict::AbstractDict, frames_at::AbstractDict)
    return Dict(ip => (i=findlast(fr -> fr.owned, frames_at[ip]); i === nothing ? frames : frames[1:i])
    for (ip, frames) in lidict)
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
    return @sprintf("%.2f", 100p), @sprintf("%.2f", 100 * sqrt(p * (1 - p) / n))
end

"`bytes` with a unit, for a table."
format_bytes(bytes::Integer) = bytes >= 2^20 ? @sprintf("%.2f MB", bytes / 2^20) : @sprintf("%.2f kB", bytes / 2^10)

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

"Write `counts.tsv` (`table`, `row`, `value`), which `compare` reads."
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
             @sprintf("%.2f", ms_per_run * get(t.onstack, k, 0) / n)] for k in rows_of(names));
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

"An optional package from the default environment (PProf, JET); not a dependency of this script."
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

"The configuration and fields every mode runs. One year is run first so compilation is not profiled; `greb_model!` restores `fields`, so they are reused."
function workload(jld2_dir::AbstractString)
    isdir(jld2_dir) || error("JLD2 data directory not found: $jld2_dir. Set GREB_DATA or pass a path.")
    cfg = preset(:full_model; corrections=Stored())
    fields = load_climatology(jld2_dir; dataset=:ncep)
    without_logs() do
        greb_model!(RunSpec(ctrl=1, scnr=0), cfg; jld2_dir=jld2_dir, fields=fields)
    end
    GC.gc()
    return cfg, fields
end

"Sample `years`-year runs until `target` samples or the mode's run limit; `mode` is `:step` or `:setup`. Returns the report directory."
function profile_cpu(jld2_dir::AbstractString, mode::Symbol; years::Int=CPU_MODES[mode].years,
        target::Int=TARGET_SAMPLES, out_dir=nothing, viewer=nothing)
    years >= 1 || throw(ArgumentError("years must be at least 1, got $years"))
    viewer in (nothing, "pprof") || throw(ArgumentError("unknown viewer $(repr(viewer)); expected pprof"))
    PProf = viewer === nothing ? nothing : load_optional("PProf", "e4faabce-9ead-11e9-39d9-4379958e3056")   # before the run, so a missing package fails early
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
        elapsed += @elapsed without_logs() do
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

    dir, info = open_report(string(mode), WORKLOAD, out_dir)
    push!(info, "threads" => "$(Threads.nthreads()) (sampled: $(join(threads_seen, ", ")))",
        "runs" => "$runs of $years simulated year(s)",
        "elapsed" => "$(round(elapsed, digits=2)) s ($(round(ms_per_run, digits=1)) ms per run)",
        "samples" => "$n ($(round(1000 * elapsed / n, digits=2)) ms per sample)")

    frames_at = Dict(ip => Frame.(frames) for (ip, frames) in lidict)
    trimmed = trim_to_package(data, samples, frames_at)
    trimmed_lidict = trim_lidict(lidict, frames_at)
    stacks = distinct_stacks(data, samples, frames_at)
    t = tally(stacks)
    mincount = max(2, n ÷ 1000)   # rows below 0.1 percent are noise
    write_counts(joinpath(dir, "counts.tsv"), t, last(info[2]))
    # Sorted by self time, largest last.
    write_report(joinpath(dir, "flat.txt"), trimmed, trimmed_lidict; format=:flat, sortedby=:overhead, mincount)
    write_report(joinpath(dir, "tree.txt"), trimmed, trimmed_lidict; format=:tree, mincount, noisefloor=2)
    write_folded(joinpath(dir, "stacks.folded"), stacks)
    PProf === nothing || Base.invokelatest(PProf.pprof, trimmed, trimmed_lidict; web=false,
        out=joinpath(dir, "profile.pb.gz"), full_signatures=false)

    close_report(dir, "Profile: $mode", info) do io, top
        write_sampled(io, t; mincount, ms_per_run, top)
    end
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

"Record every allocation of one `years`-year run; returns the report directory."
function profile_allocs(jld2_dir::AbstractString; years::Int=1, out_dir=nothing)
    years >= 1 || throw(ArgumentError("years must be at least 1, got $years"))
    cfg, fields = workload(jld2_dir)

    Profile.Allocs.clear()
    without_logs() do
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
    amount(row) = string(row.count, ", ", format_bytes(row.bytes))
    largest_first(rows) = sort!(collect(rows); by=kv -> (-last(kv).bytes, first(kv)))
    pct(row) = @sprintf("%.2f", 100 * row.bytes / total.bytes)

    dir, info = open_report("allocs", WORKLOAD, out_dir)
    push!(info, "runs" => "1 of $years simulated year(s)",
        "allocations" => amount(total),
        "in loop" => amount(in_loop) * " (with `time_loop!` on the stack)")

    close_report(dir, "Allocations", info) do io, top
        rows_of(items) = top === nothing ? items : first(items, top)
        println(io, "## By function\n")
        println(io, "Every allocation of one run, charged to the innermost package frame of its stack.\n")
        md_table(io, ["Function", "Bytes", "Share %", "Count"],
            ([code(name), format_bytes(row.bytes), pct(row), string(row.count)] for (name, row) in rows_of(largest_first(by_func))))
        println(io, "## By line\n")
        println(io, "The ", top === nothing ? 40 : top, " largest, with the type holding the most bytes there.\n")
        md_table(io, ["Line", "Bytes", "Share %", "Count", "Largest type"],
            ([code(line), format_bytes(row.bytes), pct(row), string(row.count), first(largest(row.types, 1), 80)]
            for (line, row) in first(largest_first(by_line), top === nothing ? 40 : top)); left=(1, 5))
    end
    return dir
end

"Check one model step for runtime dispatch with JET (static, nothing is timed); returns the report directory."
function profile_dispatch(jld2_dir::AbstractString; out_dir=nothing)
    JET = load_optional("JET", "c3a54625-cd67-489e-a8e7-0a5a0ff4e31b")   # before the set-up, so a missing package fails early
    isdir(jld2_dir) || error("JLD2 data directory not found: $jld2_dir. Set GREB_DATA or pass a path.")
    (; fields, r, CO2, state, ws, Ts, Ta, To, q) = stage_inputs(jld2_dir, preset(:full_model))
    acc = GREBClimate.MonthlyAccumulator()
    records = GREBClimate.MonthlyRecord[]
    times = GREBClimate.RecordTime[]
    time = ModelTime(GREBClimate.control, 1970, 1)

    step = () -> GREBClimate.time_loop!(time, CO2, Ts, Ta, To, q, records, times, fields, state, ws, acc, r)
    result = Base.invokelatest(JET.report_opt, step, (); target_modules=(GREBClimate,))
    n = length(Base.invokelatest(JET.get_reports, result))

    dir, info = open_report("dispatch", "one time_loop! step of full_model", out_dir)
    push!(info, "reports" => string(n))
    close_report(dir, "Runtime dispatch", info) do io, _
        println(io, "JET optimization analysis of one `time_loop!` step, reports from the package only. ",
            "A report is a call that is dispatched at run time; none is the expected state.\n")
        if n > 0
            println(io, "```")
            Base.invokelatest(show, IOContext(io, :color => false), MIME"text/plain"(), result)
            println(io, "\n```")
        end
    end
    return dir
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

"Write one table of `compare`: each row's share in runs `a` and `b`, the difference and its standard error."
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

# ----- Command line -----

const _PROFILE_MODES = ("step", "setup", "allocs", "dispatch", "compare")

# The flags each mode takes.
const _MODE_FLAGS = (step=("years", "samples", "out", "viewer"), setup=("years", "samples", "out", "viewer"),
    allocs=("years", "out"), dispatch=("out",), compare=())

"Run the profile named by `args`; returns the exit code."
function main(args::Vector{String})
    mode, rest = pick_mode(args, _PROFILE_MODES, "step")
    flags, rest = split_flags(rest)
    check_flags(mode, flags, _MODE_FLAGS[Symbol(mode)])
    if mode == "compare"
        length(rest) == 2 || error("compare takes two report directories: the run before and the run after")
        compare_runs(rest[1], rest[2])
        return 0
    end
    length(rest) <= 1 || error("the $mode mode takes [jld2_dir], got: ", join(rest, " "))

    jld2_dir = !isempty(rest) ? rest[1] : default_data_dir()
    out_dir = get(flags, "out", nothing)
    years = haskey(flags, "years") ? parse_nonneg_int("years", flags["years"]; min=1) : nothing
    if mode == "allocs"
        profile_allocs(jld2_dir; years=something(years, 1), out_dir)
    elseif mode == "dispatch"
        profile_dispatch(jld2_dir; out_dir)
    else
        target = haskey(flags, "samples") ? parse_nonneg_int("samples", flags["samples"]; min=1) : TARGET_SAMPLES
        profile_cpu(jld2_dir, Symbol(mode); years=something(years, CPU_MODES[Symbol(mode)].years), target, out_dir,
            viewer=get(flags, "viewer", nothing))
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
