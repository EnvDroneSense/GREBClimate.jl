# Which test files exist, which CI shard runs each, and how a developer picks a
# subset. Used by runtests.jl and by live.jl.
#
# `SHARD` is (file, CI shard, smallest tier that includes it, cost). Shards are
# grouped by measured runtime: heavy is the greb_model! integration suite, the
# golden regression and Aqua; light is the rest. Tiers are cumulative:
#   smoke    cheap files, about 40 s
#   standard + the kernel, budget and preset files
#   full     everything (the pre-commit and CI gate)
# `cost` is the approximate seconds of test time under Pkg.test. test/parallel.jl
# uses it only to balance the groups it runs side by side, so a stale value
# costs balance, not correctness.
const SHARD = [
    ("test_config.jl",     "light", "smoke",      8),
    ("test_calendar.jl",   "light", "smoke",      3),
    ("test_processes.jl",  "light", "smoke",      3),
    ("test_scenario.jl",   "light", "smoke",      3),
    ("test_state.jl",      "light", "smoke",      3),
    ("test_output.jl",     "light", "smoke",     10),
    ("test_io.jl",         "light", "smoke",      9),
    ("test_invariants.jl", "light", "smoke",     10),
    ("test_presets.jl",    "light", "standard",  16),
    ("test_budgets.jl",    "light", "standard",  11),
    ("test_experiment_tools.jl", "light", "standard", 5),
    ("test_physics.jl",    "light", "standard",  18),
    ("test_threading.jl",  "light", "full",      25),
    ("test_model.jl",      "heavy", "full",      69),
    ("test_ensemble.jl",   "heavy", "full",      30),
    ("test_golden.jl",     "heavy", "full",      21),
    ("test_aqua.jl",       "heavy", "full",      13),
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
