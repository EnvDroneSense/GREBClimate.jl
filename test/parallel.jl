# Run the test files in several Julia processes at once, so the full suite takes
# the time of its slowest group instead of the sum of all files.
#
#   julia --project=. test/parallel.jl [--jobs=4] [--tier=full] [--shard=light|heavy]
#
# Each group is one `Pkg.test(test_args = files)`, so it runs with
# --check-bounds=yes like the serial gate. Groups are balanced by the `cost`
# column of `SHARD`. test_threading.jl times subprocesses, so it runs alone after
# the groups finish. Output goes to one log per group; the logs of failed groups
# are tailed here. Exit status is non-zero if any group failed.
#
# Each process compiles and holds its own model fields (about 1 GB each), and
# the cores are shared, so do not set --jobs above the core count divided by
# three.

include(joinpath(@__DIR__, "support", "selection.jl"))

const ROOT = normpath(joinpath(@__DIR__, ".."))
const SERIAL = ("test_threading.jl",)

function parse_options(args)
    opts = Dict("jobs" => "4", "tier" => "full", "shard" => "all")
    for a in args
        m = match(r"^--(jobs|tier|shard)=(.+)$", a)
        m === nothing && error("unknown argument `$a`; use --jobs=N, --tier=NAME, --shard=light|heavy")
        opts[m.captures[1]] = m.captures[2]
    end
    jobs = tryparse(Int, opts["jobs"])
    (jobs === nothing || jobs < 1) && error("--jobs must be a positive integer, got `$(opts["jobs"])`")
    return jobs, opts["tier"], opts["shard"]
end

"Greedy longest-first packing of `files` into `n` groups of near-equal cost."
function balance(files, n)
    cost = Dict(row[1] => row[4] for row in SHARD)
    groups = [String[] for _ in 1:n]
    load = zeros(n)
    for f in sort(files; by = f -> -cost[f])
        i = argmin(load)
        push!(groups[i], f)
        load[i] += cost[f]
    end
    return filter!(!isempty, groups), cost
end

function launch(files, log)
    code = "using Pkg; Pkg.test(test_args = $(repr(files)))"
    cmd = `$(Base.julia_cmd()) --project=$ROOT -e $code`
    return run(pipeline(ignorestatus(cmd); stdout = log, stderr = log); wait = false)
end

function summary_line(log)
    lines = readlines(log)
    i = findlast(l -> startswith(l, "GREBClimate.jl"), lines)
    return i === nothing ? "no summary line" : strip(lines[i])
end

function main(args)
    jobs, tier, shard = parse_options(args)
    selected = select_tests(; tier = tier == "full" ? "" : tier, shard = shard)
    parallel = filter(!in(SERIAL), selected)
    serial = filter(in(SERIAL), selected)
    groups, cost = balance(parallel, min(jobs, length(parallel)))
    for s in serial
        push!(groups, [s])
    end

    logdir = mktempdir(; cleanup = false)
    failed = String[]
    started = time()
    println("running ", length(selected), " files in ", length(groups) - length(serial),
            " parallel groups", isempty(serial) ? "" : " then threading alone", "; logs in ", logdir)

    results = NamedTuple[]
    function finish(i, group, proc, t0)
        wait(proc)
        log = joinpath(logdir, "group$i.log")
        ok = success(proc)
        ok || push!(failed, log)
        push!(results, (group = i, files = join(replace.(group, r"^test_|\.jl$" => ""), " "),
                        seconds = round(time() - t0; digits = 1), ok = ok, summary = summary_line(log)))
    end

    n_par = length(groups) - length(serial)
    procs = [(i, g, launch(g, joinpath(logdir, "group$i.log")), time()) for (i, g) in enumerate(groups[1:n_par])]
    # One task per process, so each group's time stops when that group does
    @sync for p in procs
        @async finish(p...)
    end
    for (k, g) in enumerate(groups[n_par+1:end])
        i = n_par + k
        finish(i, g, launch(g, joinpath(logdir, "group$i.log")), time())
    end

    for r in sort(results; by = r -> r.group)
        println("group ", r.group, "  ", r.ok ? "ok    " : "FAILED", "  ", lpad(r.seconds, 6), " s  ",
                r.files, "\n          ", r.summary)
    end
    println("wall ", round(time() - started; digits = 1), " s; this proves only the files listed above")
    for log in failed
        println("\n--- tail of ", log)
        foreach(println, last(readlines(log), 40))
    end
    exit(isempty(failed) ? 0 : 1)
end

main(ARGS)
