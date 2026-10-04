# Package hygiene: stale deps, missing compat entries, method ambiguities,
# type piracy, unbound type parameters. In the heavy shard: Aqua costs ~24 s.
# persistent_tasks is skipped: __init__ only registers a DataDep and starts
# nothing, and the check costs 21 s.

using Aqua

@testset "Aqua quality assurance" begin
    Aqua.test_all(GREBClimate; persistent_tasks = false)
end
