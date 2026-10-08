# Which test files exist, which CI shard runs each, and how a developer picks a
# subset. Used by runtests.jl and by live.jl.
#
# `SHARD` is (file, CI shard, smallest tier that includes it). Shards are
# grouped by measured runtime: heavy is the greb_model! integration suite, the
# golden regression and Aqua; light is the rest. Tiers are cumulative:
#   smoke    cheap files, about 40 s
#   standard + the kernel, budget and preset files
#   full     everything (the pre-commit and CI gate)
const SHARD = [
    ("test_config.jl",     "light", "smoke"),
    ("test_calendar.jl",   "light", "smoke"),
    ("test_processes.jl",  "light", "smoke"),
    ("test_scenario.jl",   "light", "smoke"),
    ("test_state.jl",      "light", "smoke"),
    ("test_output.jl",     "light", "smoke"),
    ("test_io.jl",         "light", "smoke"),
    ("test_invariants.jl", "light", "smoke"),
    ("test_presets.jl",    "light", "standard"),
    ("test_budgets.jl",    "light", "standard"),
    ("test_experiment_tools.jl", "light", "standard"),
    ("test_physics.jl",    "light", "standard"),
    ("test_threading.jl",  "light", "full"),
    ("test_model.jl",      "heavy", "full"),
    ("test_ensemble.jl",   "heavy", "full"),
    ("test_golden.jl",     "heavy", "full"),
    ("test_aqua.jl",       "heavy", "full"),
]

const TIERS = ["smoke", "standard", "full"]

"""
`"state"`, `"test_state"` and `"test_state.jl"` all name test_state.jl. Throws
on a name that is not in `SHARD`, so a typo cannot silently run nothing.
"""
function resolve_test_file(name)
    base = replace(strip(name), r"\.jl$" => "")
    startswith(base, "test_") || (base = "test_" * base)
    file = base * ".jl"
    any(row -> row[1] == file, SHARD) ||
        error("no test file `$name` in SHARD; known: ", join((r[1] for r in SHARD), ", "))
    return file
end

"""
The files to run, in `SHARD` order. `names` (from `GREB_TEST_FILES` or the
`Pkg.test(test_args = ...)` arguments) wins over `tier`, which wins over
`shard`. All empty means everything.
"""
function select_tests(; names = String[], tier = "", shard = "all")
    if !isempty(names)
        wanted = Set(resolve_test_file.(names))
        return [row[1] for row in SHARD if row[1] in wanted]
    elseif !isempty(tier)
        rank = findfirst(==(tier), TIERS)
        rank === nothing && error("GREB_TEST_TIER must be one of ", join(TIERS, ", "), ", got `$tier`")
        return [row[1] for row in SHARD if findfirst(==(row[3]), TIERS) <= rank]
    else
        return [row[1] for row in SHARD if shard == "all" || shard == row[2]]
    end
end
