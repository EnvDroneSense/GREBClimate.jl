# run_ensemble: members side by side give the results of running them alone.

@testset "run_ensemble: every member equals its own sequential run" begin
    fields = synthetic_fields()
    untouched = deepcopy(fields)
    run = RunSpec(ctrl = 1, scnr = 1)
    configs = [preset(:co2_double; corrections = NoCorrections()),
               preset(:co2_double; corrections = NoCorrections(), hydrology = (evaporation = :skin,)),
               preset(:co2_double; corrections = SpinUp(1))]
    alone = [run_synthetic(run, c; fields = deepcopy(fields)) for c in configs]
    same(a, b) = length(a) == length(b) &&
                 all(isequal(getfield(x, f), getfield(y, f)) for (x, y) in zip(a, b) for f in fieldnames(typeof(x)))

    # Fewer copies of the fields than members, so a copy is reused; no log line gets out
    together = @test_logs run_ensemble(run, configs; fields, ntasks = 2, allow_uninitialized = true)
    @test length(together) == length(configs)
    for (a, b) in zip(alone, together)
        @test same(a.ctrl, b.ctrl) && same(a.scnr, b.scnr)
    end
    # The members differ from each other, so the comparison above means something
    @test !same(alone[1].scnr, alone[2].scnr) && !same(alone[1].ctrl, alone[3].ctrl)

    # The caller's fields are not changed
    @test all(isequal(getfield(fields, f), getfield(untouched, f)) for f in fieldnames(ClimateFields))
end

@testset "run_ensemble: the copies of the fields are released when it returns" begin
    fields = synthetic_fields()
    run = RunSpec(ctrl = 0, scnr = 1)
    config = preset(:co2_double; corrections = NoCorrections())
    live() = (GC.gc(true); GC.gc(true); Base.gc_live_bytes())
    quiet() do
        greb_model!(run, config; fields = deepcopy(fields), allow_uninitialized = true)   # compile first
    end
    before = live()
    run_ensemble(r -> length(r.scnr), run, fill(config, 3); fields, ntasks = 3, allow_uninitialized = true)
    @test live() - before < Base.summarysize(fields) ÷ 2   # three copies were made, none is kept
end

@testset "run_ensemble: logging, errors and arguments" begin
    fields = synthetic_fields()
    run = RunSpec(ctrl = 1, scnr = 0)
    configs = [preset(:full_model; corrections = NoCorrections()) for _ in 1:3]
    ensemble(f = identity; kwargs...) = run_ensemble(f, run, configs; fields, allow_uninitialized = true, kwargs...)

    # With a logger the members' lines come through. `reduce` runs per member;
    # more tasks than members is the same as one each
    means = @test_logs (:info, r"Control run") match_mode = :any ensemble(
        result -> global_mean(result.ctrl[end].Ts); logger = Base.current_logger(), ntasks = 50)
    @test means isa Vector{Float64} && length(means) == 3 && allequal(means) && isfinite(means[1])

    # A failing member is reported after the others have finished, and nothing is left blocked
    finished = Threads.Atomic{Int}(0)
    fails_once = let first = Threads.Atomic{Bool}(true)
        result -> Threads.atomic_xchg!(first, false) ? error("member failed") : Threads.atomic_add!(finished, 1)
    end
    @test_throws CompositeException ensemble(fails_once; ntasks = 2)
    @test finished[] == 2

    @test_throws ArgumentError ensemble(ntasks = 0)
    @test isempty(run_ensemble(run, Config[]; fields))
end
