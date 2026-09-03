# Package hygiene: stale deps, missing compat entries, method ambiguities,
# type piracy, unbound type parameters. In the heavy shard because Aqua costs
# ~24s (ambiguities 14s, stale_deps 8s, the rest under 1s).
#
# persistent_tasks is skipped: it re-precompiles the package in a subprocess to
# find tasks or timers left running, and __init__ only registers a DataDep - it
# starts nothing. Not worth 21s.

using Aqua

@testset "Aqua quality assurance" begin
    Aqua.test_all(GREBClimate; persistent_tasks = false)
end
