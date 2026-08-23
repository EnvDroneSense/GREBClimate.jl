# Threaded-vs-serial equivalence (spawns its own subprocesses).

@testset "threaded circulation matches serial (subprocess -t 1 vs -t 2)" begin
    # `tendencies!` runs circulation!(Ta) and circulation!(q) concurrently only
    # when `Threads.nthreads() > 1` AND `ws_a !== ws_q` (see src/tendencies.jl).
    # Thread count is fixed at Julia startup, so a single-threaded `Pkg.test()`
    # can never reach that branch - it went untested until 2026-08-21. Spawning
    # both counts explicitly keeps this honest however the suite is invoked.
    utils = joinpath(@__DIR__, "testutils.jl")
    script = """
        using GREBClimate
        using Test
        include(raw"$(utils)")
        cfg = create_experiment_config(:full_model)
        result = quiet() do
            greb_model!(RunSpec(flux = 0, ctrl = 1, scnr = 0), cfg;
                        jld2_dir = "", fields = synthetic_fields(),
                        allow_uninitialized = true)
        end
        print(Threads.nthreads())
        for rec in result.ctrl
            print(" ", gmean(rec.Ts), " ", gmean(rec.Ta), " ", gmean(rec.q))
        end
    """
    # Only the executable from julia_cmd(), not its flags: under Pkg.test those
    # include --check-bounds=yes, which would force the subprocess to recompile
    # the world and make this test ~10x slower.
    exe = first(Base.julia_cmd())
    project = normpath(joinpath(@__DIR__, ".."))
    run_at(n) = begin
        cmd = `$exe --startup-file=no --project=$project -t $n -e $script`
        parts = split(strip(read(cmd, String)))
        (nthreads = parse(Int, parts[1]), digest = parse.(Float64, parts[2:end]))
    end

    serial = run_at(1)
    threaded = run_at(2)

    # the subprocesses really did run at the requested thread counts
    @test serial.nthreads == 1
    @test threaded.nthreads == 2
    # 12 months x 3 quantities
    @test length(serial.digest) == 36
    @test length(threaded.digest) == length(serial.digest)

    @test threaded.digest == serial.digest
end
