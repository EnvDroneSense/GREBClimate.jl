using GREBClimate
using Test

include(joinpath("support", "testutils.jl"))

# One file per subject; `SHARD` decides which CI job runs each. Set
# GREB_TEST_SHARD=light|heavy to run one group; unset (or "all") runs
# everything. Shards are grouped by measured runtime: heavy is the greb_model!
# integration suite, the golden regression and Aqua; light is the rest.
const SHARD = [
    ("test_config.jl",     "light"),
    ("test_presets.jl",    "light"),
    ("test_processes.jl",  "light"),
    ("test_scenario.jl",   "light"),
    ("test_state.jl",      "light"),
    ("test_output.jl",     "light"),
    ("test_budgets.jl",    "light"),
    ("test_physics.jl",    "light"),
    ("test_io.jl",         "light"),
    ("test_invariants.jl", "light"),
    ("test_threading.jl",  "light"),
    ("test_model.jl",      "heavy"),
    ("test_golden.jl",     "heavy"),
    ("test_aqua.jl",       "heavy"),
]

shard = get(ENV, "GREB_TEST_SHARD", "all")
@testset "GREBClimate.jl" begin
    for (file, group) in SHARD
        (shard == "all" || shard == group) || continue
        include(file)
    end
end
