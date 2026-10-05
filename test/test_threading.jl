# Threaded-vs-serial equivalence (spawns its own subprocesses).

@testset "threaded circulation matches serial (subprocess -t 1 vs -t 2)" begin
    # `tendencies!` runs the two circulation! calls concurrently only when
    # `Threads.nthreads() > 1` and `ws_a !== ws_q`. Thread count is fixed at
    # Julia startup, so both counts are started as subprocesses.
    utils = joinpath(@__DIR__, "support", "testutils.jl")
    script = """
        using GREBClimate
        using Test
        include(raw"$(utils)")
        cfg = preset(:full_model; corrections = NoCorrections())
        result = quiet() do
            greb_model!(RunSpec(ctrl = 1, scnr = 0), cfg;
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

    @test isequal(threaded.digest, serial.digest)
end

@testset "threaded run matches serial on the real dataset (subprocess -t 1 vs -t 2)" begin
    # The synthetic run above uses made-up fields with no sunlight, so it
    # compares little real physics. This one runs a flux-correction spin-up
    # and a control year on the dataset and compares every record field.
    if !isdir(DATA_DIR)
        @test_skip "greb_input_data/ not present"
    else
        script = """
            using GREBClimate
            result = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                fields = load_climatology(raw"$(DATA_DIR)"; dataset = :ncep)
                greb_model!(RunSpec(ctrl = 1, scnr = 0), preset(:full_model; corrections = SpinUp(1));
                            jld2_dir = raw"$(DATA_DIR)", fields = fields)
            end
            finite = all(r -> all(isfinite, r.Ts) && all(isfinite, r.Ta), result.ctrl)
            print(Threads.nthreads(), " ", finite, " ",
                  hash([getfield(r, v) for r in result.ctrl for v in fieldnames(MonthlyRecord)]))
        """
        exe = first(Base.julia_cmd())
        project = normpath(joinpath(@__DIR__, ".."))
        run_at(n) = split(strip(read(`$exe --startup-file=no --project=$project -t $n -e $script`, String)))

        serial = run_at(1)
        threaded = run_at(2)
        @test serial[1:2] == ["1", "true"]
        @test threaded[1:2] == ["2", "true"]
        @test threaded[3] == serial[3]
    end
end
