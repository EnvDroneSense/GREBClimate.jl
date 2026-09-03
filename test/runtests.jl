using GREBClimate
using Test

include("testutils.jl")

# One file per subject area; `SHARD` decides which CI job runs each. Set
# GREB_TEST_SHARD=light|heavy to run one group, or leave it unset (or "all")
# so plain `Pkg.test()`/`]test` runs everything - sharding is opt-in, CI-only.
#
# Shards are balanced by measured runtime, not by file count: the heavy job is
# the greb_model! integration suite plus the golden regression, the light job is
# everything else plus the threading subprocesses.
const SHARD = [
    ("test_config.jl",     "light"),
    ("test_state.jl",      "light"),
    ("test_output.jl",     "light"),
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
