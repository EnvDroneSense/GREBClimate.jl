# Sampling profile of GREBClimate.jl: where the time of a run goes.
#
#   julia --project=. -t 1 benchmark/profile.jl [mode] [jld2_dir] [--years=N] [--out=DIR]
#
#   step - profile a control run of `--years` simulated years (default 50)
#
# Writes header.txt, flat.txt and tree.txt to
# benchmark/profiles/<date>-<commit>-<mode>/ (or `--out`). Run it at `-t 1`: on
# Windows the sampler records the first thread only, so work on other threads
# shows up as waiting.

using GREBClimate
using Dates
using Profile

include("common.jl")

const MIN_SAMPLES = 5_000        # below this only the largest rows mean anything
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
Profile a `years`-year `:full_model` control run on the stored flux corrections
(no spin-up) and write the reports to `out_dir`; returns `out_dir`.
"""
function profile_step(jld2_dir::AbstractString; years::Int=50, out_dir=nothing)
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

    fields_r = deepcopy(fields)
    Profile.clear()
    Profile.init(n=PROFILE_BUFFER, delay=0.001)
    elapsed = @elapsed quiet() do
        @profile greb_model!(RunSpec(ctrl=years, scnr=0), cfg; jld2_dir=jld2_dir, fields=fields_r)
    end
    Profile.is_buffer_full() && error("the profile buffer filled up; lower --years")

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
        years:       $years
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
    out_dir = get(flags, "out", nothing)

    if mode == "step"
        profile_step(jld2_dir; years, out_dir)
    end
end
