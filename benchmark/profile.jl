# Sampling profile of GREBClimate.jl: where the time of a run goes.
#
#   julia --project=. -t 1 benchmark/profile.jl [mode] [jld2_dir] [--years=N] [--samples=N] [--out=DIR]
#
#   step - profile control runs of `--years` simulated years (default 50), repeated
#          until `--samples` samples are collected (default 10000)
#
# Writes header.txt, category.txt, owned.txt, flat.txt and tree.txt to
# benchmark/profiles/<date>-<commit>-<mode>/ (or `--out`). Run it at `-t 1`: on
# Windows the sampler records the first thread only, so work on other threads
# shows up as waiting.

using GREBClimate
using Dates
using Printf
using Profile

include("common.jl")

const MIN_SAMPLES = 5_000        # below this only the largest rows mean anything
const MAX_RUNS = 8               # the sampler's rate varies, so a run is repeated, up to this often
const PROFILE_BUFFER = 10_000_000
const META_WORDS = 6             # threadid, taskid, cycle clock, sleep state, 0, 0 close each sample

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

"True if any frame at instruction pointer `ip` lies in the package's `src/`."
function in_package(ip::UInt64, lidict, srcdir::AbstractString)
    any(fr -> startswith(normpath(string(fr.file)), srcdir), lidict[ip])
end

"""
Copy of `data` with the frames outside the package's outermost frame removed
from every sample, so a tree starts at the model and not at Julia's start-up.
A sample with no package frame is kept whole.
"""
function trim_to_package(data::Vector{UInt64}, lidict, samples)
    srcdir = normpath(joinpath(REPO, "src"))
    owned = Dict{UInt64,Bool}()
    trimmed = UInt64[]
    sizehint!(trimmed, length(data))
    for (frames, meta) in samples
        # Frames are stored innermost first, so the outermost package frame is the last match.
        last_owned = findlast(i -> get!(() -> in_package(data[i], lidict, srcdir), owned, data[i]), frames)
        keep = last_owned === nothing ? frames : first(frames):frames[last_owned]
        append!(trimmed, view(data, keep))
        append!(trimmed, view(data, meta))
    end
    return trimmed
end

"One stack frame, reduced to what the reports need."
struct Frame
    func::String
    file::String
    line::Int
    from_c::Bool
    owned::Bool      # lies in the package's `src/`
end

"`name` without the `#name#12` wrapping of keyword bodies and closures, whose number changes between builds."
function clean_name(name::AbstractString)
    m = match(r"^#+([^#]+)#+\d*$", name)
    (m === nothing || all(isdigit, m.captures[1])) ? String(name) : String(m.captures[1])
end

"The stack of every sample as `Frame`s, innermost first, inlined frames included."
function sample_stacks(data::Vector{UInt64}, lidict, samples)
    srcdir = normpath(joinpath(REPO, "src"))
    cache = Dict{UInt64,Vector{Frame}}()
    return map(samples) do (frames, _)
        stack = Frame[]
        for i in frames
            append!(stack, get!(cache, data[i]) do
                [Frame(clean_name(string(fr.func)), replace(string(fr.file), '\\' => '/'), fr.line,
                     fr.from_c, startswith(normpath(string(fr.file)), srcdir)) for fr in lidict[data[i]]]
            end)
        end
        stack
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

"`count` of `n` samples as a percentage and its standard error in percentage points, as padded strings."
function share(count::Int, n::Int)
    p = count / n
    return @sprintf("%7.2f", 100p), @sprintf("%6.2f", 100 * sqrt(p * (1 - p) / n))
end

"The `top` largest entries of `counts` as `name (count)`, joined by commas."
function largest(counts::Dict{String,Int}, top::Int=4)
    rows = first(sort!(collect(counts); by=last, rev=true), top)
    return join((string(k, " (", v, ")") for (k, v) in rows), ", ")
end

"""
Write `owned.txt`: samples per package function and per package line, where a
sample belongs to its innermost package frame, so library time is charged to
the line that called it.
"""
function write_owned(path::AbstractString, stacks::Vector{Vector{Frame}}; mincount::Int)
    n = length(stacks)
    self = Dict{String,Int}()
    onstack = Dict{String,Int}()
    lines = Dict{Tuple{String,Int,String},Int}()
    for stack in stacks
        o = owner(stack)
        if o === nothing
            self["(outside the package)"] = get(self, "(outside the package)", 0) + 1
            continue
        end
        fr, func = o
        self[func] = get(self, func, 0) + 1
        key = (basename(fr.file), fr.line, func)
        lines[key] = get(lines, key, 0) + 1
        for name in unique!([f.func for f in stack if f.owned && f.func != "macro expansion"])
            onstack[name] = get(onstack, name, 0) + 1
        end
    end

    open(path, "w") do io
        println(io, "Samples: ", n, ". Shares in percent, +/- is one standard error in percentage points.")
        println(io, "Rows below ", mincount, " samples are left out of the tables but counted in the total.\n")
        println(io, "By function. self: the function holds the innermost package frame. on stack: it is anywhere in the stack.")
        println(io, "   self  share%    +/-  on stack  share%  function")
        names = sort!(collect(union(keys(self), keys(onstack))); by=k -> (get(self, k, 0), get(onstack, k, 0)), rev=true)
        for name in names
            s, t = get(self, name, 0), get(onstack, name, 0)
            max(s, t) >= mincount || continue
            p, se = share(s, n)
            println(io, lpad(s, 7), " ", p, " ", se, lpad(t, 10), " ", share(t, n)[1], "  ", name)
        end
        println(io, lpad(sum(values(self)), 7), " ", share(sum(values(self)), n)[1], "         total\n")

        println(io, "By line: the innermost package frame. A `@turbo` loop is one row, at its `@turbo for` line.")
        println(io, "   self  share%    +/-  line")
        for ((file, line, func), c) in sort!(collect(lines); by=last, rev=true)
            c >= mincount || break
            p, se = share(c, n)
            println(io, lpad(c, 7), " ", p, " ", se, "  ", file, ":", line, "  ", func)
        end
    end
end

"""
Write `category.txt`: samples per kind of cost, judged by the innermost frame,
with the largest frames and the package functions they are charged to.
"""
function write_category(path::AbstractString, stacks::Vector{Vector{Frame}})
    n = length(stacks)
    count = Dict{String,Int}()
    frames = Dict{String,Dict{String,Int}}()
    owners = Dict{String,Dict{String,Int}}()
    for stack in stacks
        isempty(stack) && continue
        fr = first(stack)
        class = classify(fr)
        count[class] = get(count, class, 0) + 1
        label = string(fr.func, " ", basename(fr.file), ":", fr.line)
        d = get!(Dict{String,Int}, frames, class)
        d[label] = get(d, label, 0) + 1
        o = owner(stack)
        name = o === nothing ? "(outside the package)" : o[2]
        d = get!(Dict{String,Int}, owners, class)
        d[name] = get(d, name, 0) + 1
    end

    classes = sort!(collect(count); by=last, rev=true)
    open(path, "w") do io
        println(io, "Samples: ", n, ". Each sample is classed by its innermost frame.")
        println(io, "Shares in percent, +/- is one standard error in percentage points.\n")
        println(io, "samples  share%    +/-  class")
        for (class, c) in classes
            p, se = share(c, n)
            println(io, lpad(c, 7), " ", p, " ", se, "  ", class)
        end
        println(io, lpad(sum(values(count)), 7), " ", share(sum(values(count)), n)[1], "         total")
        for (class, _) in classes
            println(io, "\n", class)
            println(io, "  frames:     ", largest(frames[class]))
            println(io, "  charged to: ", largest(owners[class]))
        end
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

"A directory for this run's reports under `benchmark/profiles/`, numbered if the name is taken."
function default_out_dir(mode::AbstractString, commit::AbstractString, dirty::Bool)
    name = string(Dates.today(), "-", commit, dirty ? "-dirty" : "", "-", mode)
    base = joinpath(REPO, "benchmark", "profiles", name)
    dir, n = base, 1
    while ispath(dir)
        n += 1
        dir = string(base, "-", n)
    end
    return dir
end

"Write `Profile.print` output for `data` to `path` without truncating long lines."
function write_report(path::AbstractString, data, lidict; kwargs...)
    open(path, "w") do io
        Profile.print(IOContext(io, :displaysize => (200, 250)), data, lidict; kwargs...)
    end
end

"""
Profile `years`-year `:full_model` control runs on the stored flux corrections
(no spin-up), repeated until `target` samples are collected, and write the
reports to `out_dir`; returns `out_dir`.
"""
function profile_step(jld2_dir::AbstractString; years::Int=50, target::Int=10_000, out_dir=nothing)
    isdir(jld2_dir) || error("JLD2 data directory not found: $jld2_dir. Set GREB_DATA or pass a path.")
    years >= 1 || throw(ArgumentError("years must be at least 1, got $years"))

    if Threads.nthreads() > 1
        @warn "Profiling with $(Threads.nthreads()) threads: on Windows only the first thread is sampled. Use -t 1."
    end

    cfg = preset(:full_model; corrections=Stored())
    fields = load_greb_jld2!(jld2_dir; dataset=:ncep)
    # The run's progress lines are not part of the workload.
    quiet(f) = redirect_stdout(devnull) do
        Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
    end

    # Warm-up, so compilation is not sampled.
    quiet() do
        greb_model!(RunSpec(ctrl=1, scnr=0), cfg; jld2_dir=jld2_dir, fields=deepcopy(fields))
    end

    Profile.clear()
    Profile.init(n=PROFILE_BUFFER, delay=0.001)
    elapsed = 0.0
    runs = 0
    while runs < MAX_RUNS && (runs == 0 || length(split_samples(Profile.fetch(include_meta=true))) < target)
        fields_r = deepcopy(fields)  # a fresh copy per run, outside the profiled region
        elapsed += @elapsed quiet() do
            @profile greb_model!(RunSpec(ctrl=years, scnr=0), cfg; jld2_dir=jld2_dir, fields=fields_r)
        end
        runs += 1
        Profile.is_buffer_full() && error("the profile buffer filled up; lower --years or --samples")
    end

    data, lidict = Profile.retrieve(include_meta=true)
    samples = split_samples(data)
    nsamples = length(samples)
    nsamples > 0 || error("the profiler recorded no samples")
    threads_seen = sort!(unique(Int(data[first(meta)]) for (_, meta) in samples))

    commit, dirty = repo_state()
    dir = something(out_dir, default_out_dir("step", commit, dirty))
    mkpath(dir)

    header = """
        mode:        step (full_model control, stored corrections)
        commit:      $commit$(dirty ? " (tracked files modified)" : "")
        date:        $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM"))
        julia:       $VERSION
        threads:     $(Threads.nthreads()) (sampled: $(join(threads_seen, ", ")))
        years:       $(runs * years) ($runs runs of $years)
        elapsed:     $(round(elapsed, digits=2)) s
        samples:     $nsamples
        interval:    $(round(1000 * elapsed / nsamples, digits=2)) ms per sample
        """
    write(joinpath(dir, "header.txt"), header)

    trimmed = trim_to_package(data, lidict, samples)
    mincount = max(2, nsamples ÷ 1000)   # rows below 0.1 percent are noise
    # Sorted by self time, largest last.
    write_report(joinpath(dir, "flat.txt"), trimmed, lidict; format=:flat, sortedby=:overhead, mincount)
    write_report(joinpath(dir, "tree.txt"), trimmed, lidict; format=:tree, mincount, noisefloor=2)
    stacks = sample_stacks(data, lidict, samples)
    write_owned(joinpath(dir, "owned.txt"), stacks; mincount)
    write_category(joinpath(dir, "category.txt"), stacks)

    print(header)
    nsamples >= MIN_SAMPLES ||
        @warn "Only $nsamples samples (fewer than $MIN_SAMPLES): shares of small rows are not reliable. Raise --years."
    println("reports written to ", dir)
    return dir
end

const _PROFILE_MODES = ("step",)

if abspath(PROGRAM_FILE) == @__FILE__
    mode, rest = if isempty(ARGS) || startswith(ARGS[1], "--")
        ("step", copy(ARGS))
    elseif ARGS[1] in _PROFILE_MODES
        (ARGS[1], ARGS[2:end])
    else
        error("unknown mode $(repr(ARGS[1])); expected one of $(join(_PROFILE_MODES, ", "))")
    end

    flags, rest = split_flags(rest)
    jld2_dir = !isempty(rest) ? rest[1] : default_data_dir()
    years = haskey(flags, "years") ? parse_nonneg_int("years", flags["years"]) : 50
    target = haskey(flags, "samples") ? parse_nonneg_int("samples", flags["samples"]) : 10_000
    out_dir = get(flags, "out", nothing)

    if mode == "step"
        profile_step(jld2_dir; years, target, out_dir)
    end
end
