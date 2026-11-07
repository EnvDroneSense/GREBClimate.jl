# Which test files exist, which CI shard runs each, and how a developer picks a
# subset. Used by runtests.jl and by live.jl.
#
# `SHARD` is (file, CI shard, smallest tier that includes it, cost). Shards are
# grouped by measured runtime: heavy is the greb_model! integration suite, the
# golden regression and Aqua; light is the rest. Tiers are cumulative:
#   smoke    cheap files, about 40 s
#   standard + the kernel, budget and preset files
#   full     everything (the pre-commit and CI gate)
# `cost` is the approximate seconds of test time under Pkg.test. test/support/parallel.jl
# uses it only to balance the groups it runs side by side, so a stale value
# costs balance, not correctness.
const SHARD = [
    ("test_config.jl", "light", "smoke", 8),
    ("test_calendar.jl", "light", "smoke", 3),
    ("test_processes.jl", "light", "smoke", 3),
    ("test_scenario.jl", "light", "smoke", 4),
    ("test_state.jl", "light", "smoke", 3),
    ("test_output.jl", "light", "smoke", 9),
    ("test_io.jl", "light", "smoke", 14),
    ("test_selection.jl", "light", "smoke", 1),
    ("test_invariants.jl", "light", "smoke", 10),
    ("test_presets.jl", "light", "standard", 16),
    ("test_budgets.jl", "light", "standard", 11),
    ("test_experiment_tools.jl", "light", "standard", 5),
    ("test_physics.jl", "light", "standard", 26),
    ("test_threading.jl", "light", "full", 25),
    ("test_enzyme.jl", "light", "full", 45),
    ("test_model.jl", "heavy", "full", 49),
    ("test_ensemble.jl", "heavy", "full", 30),
    ("test_golden.jl", "heavy", "full", 21),
    ("test_aqua.jl", "heavy", "full", 13)
]

# Run only by name (`Pkg.test(test_args = ["dataset_tools"])`), never as part of a
# tier, a shard or the full suite: the dataset tools change rarely, and their
# tests are needed only when they or the dataset loaders do.
const OPT_IN = ["test_dataset_tools.jl"]

const TIERS = ["smoke", "standard", "full"]

"""
`"state"`, `"test_state"` and `"test_state.jl"` all name test_state.jl. Throws
on a name that is not in `SHARD`, so a typo cannot silently run nothing.
"""
function resolve_test_file(name)
    base = replace(strip(name), r"\.jl$" => "")
    startswith(base, "test_") || (base = "test_" * base)
    file = base * ".jl"
    file in OPT_IN || any(row -> row[1] == file, SHARD) ||
        error("no test file `$name` in SHARD; known: ", join(vcat([r[1] for r in SHARD], OPT_IN), ", "))
    return file
end

"""
The files to run, in `SHARD` order. `names` (from `GREB_TEST_FILES` or the
`Pkg.test(test_args = ...)` arguments) wins over `tier`, which wins over
`shard`. All empty means everything.
"""
function select_tests(; names=String[], tier="", shard="all")
    if !isempty(names)
        wanted = Set(resolve_test_file.(names))
        return [f for f in vcat([row[1] for row in SHARD], OPT_IN) if f in wanted]
    elseif !isempty(tier)
        rank = findfirst(==(tier), TIERS)
        rank === nothing && error("GREB_TEST_TIER must be one of ", join(TIERS, ", "), ", got `$tier`")
        return [row[1] for row in SHARD if findfirst(==(row[3]), TIERS) <= rank]
    else
        return [row[1] for row in SHARD if shard == "all" || shard == row[2]]
    end
end

# Which test files a change needs: the source-to-test map of CLAUDE.md as data,
# read by test/support/changed.jl. `"@standard"` stands for every file of the standard
# tier. The first rule whose prefix starts the path wins. A path no rule names
# runs the full suite, so a new file is safe by default; test_selection.jl fails
# when a file under src/ or ext/ has no rule, so it gets one on purpose.
const SOURCE_MAP = [
    ("src/physics/", ["physics", "budgets", "invariants"]),
    ("src/core/tendencies.jl", ["physics", "budgets", "invariants"]),
    ("src/core/budgets.jl", ["budgets", "invariants"]),
    ("src/config/", ["config", "presets", "processes", "scenario"]),
    ("src/forcing/", ["config", "presets", "processes", "scenario"]),
    ("src/core/calendar.jl", ["calendar"]),
    ("src/core/state.jl", ["state", "output"]),
    ("src/core/output.jl", ["state", "output"]),
    ("src/core/postprocess.jl", ["state", "output"]),
    ("src/io.jl", ["io"]),
    ("src/data.jl", ["io"]),
    ("src/core/model.jl", ["model", "ensemble", "@standard"]),
    ("src/core/ensemble.jl", ["model", "ensemble", "@standard"]),
    ("ext/", ["enzyme"]),
    ("tools/dataset/", ["dataset_tools", "io"]),
    ("tools/experiments/", ["experiment_tools"]),
    ("tools/", ["experiment_tools", "dataset_tools"]),
]

# Paths that run the full suite although they are named: the package module,
# constants, the environments and the shared test code
const FULL_PATHS = ["src/GREBClimate.jl", "src/core/constants.jl", "Project.toml",
    "test/Project.toml", "test/support/", "test/runtests.jl", ".github/workflows/"]

# The three runner scripts in test/support/: no test file exercises them, so a
# change runs only the selection test, and the maintainer runs the script itself
const RUNNER_PATHS = ["test/support/parallel.jl", "test/support/live.jl", "test/support/changed.jl"]

# Paths no test can fail on
const NO_TEST_PATHS = ["docs/", "benchmark/", "viz/", ".claude/", "README", "test/README", "CHANGELOG", "CONTRIBUTING",
    "DATA_README", "LICENSE", ".gitignore", ".JuliaFormatter", "CITATION", "codecov"]

"""
`(files, why)` for the changed paths (relative to the repository root, with `/`):
the test files to run in `SHARD` order, and a line that says why. `why` is
`"full"` when the whole suite runs. A test file runs itself.
"""
function tests_for(paths)
    wanted = Set{String}()
    for path in paths
        if startswith(path, "test/test_") && endswith(path, ".jl")
            push!(wanted, resolve_test_file(basename(path)))
        elseif path in RUNNER_PATHS
            push!(wanted, "test_selection.jl")
        elseif any(startswith(path, p) for p in FULL_PATHS)
            return (select_tests(), "full: $path")
        elseif any(startswith(path, p) for p in NO_TEST_PATHS)
            continue
        else
            i = findfirst(r -> startswith(path, r[1]), SOURCE_MAP)
            i === nothing && return (select_tests(), "full: no rule for $path")
            for name in SOURCE_MAP[i][2]
                union!(wanted, name == "@standard" ? select_tests(tier="standard") :
                               [resolve_test_file(name)])
            end
        end
    end
    return ([f for f in vcat([row[1] for row in SHARD], OPT_IN) if f in wanted], "mapped")
end
