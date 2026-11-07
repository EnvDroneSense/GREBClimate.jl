# The source-to-test map (support/selection.jl) that test/support/changed.jl reads.

@testset "every source file has a rule; every rule names real test files" begin
    root = normpath(joinpath(@__DIR__, ".."))
    sources = [replace(relpath(joinpath(dir, f), root), "\\" => "/")
               for top in ("src", "ext") for (dir, _, fs) in walkdir(joinpath(root, top)) for f in fs if endswith(f, ".jl")]
    @test !isempty(sources)
    named(path) = any(startswith(path, p) for p in FULL_PATHS) || any(r -> startswith(path, r[1]), SOURCE_MAP)
    @test filter(!named, sources) == String[]
    for (_, names) in SOURCE_MAP, name in names
        name == "@standard" || @test resolve_test_file(name) in vcat(first.(SHARD), OPT_IN)
    end
    # Every path prefix names something that exists
    for prefix in [first.(SOURCE_MAP); FULL_PATHS[1:5]; RUNNER_PATHS]
        @test ispath(joinpath(root, prefix))
    end
end

@testset "tests_for: mapped, full and no-test changes" begin
    @test tests_for(["test/test_state.jl"]) == (["test_state.jl"], "mapped")
    @test tests_for(["src/core/calendar.jl"]) == (["test_calendar.jl"], "mapped")
    @test tests_for(["docs/src/index.md", "CHANGELOG.md"]) == (String[], "mapped")
    # The model files run their own files and the standard tier
    model, _ = tests_for(["src/core/model.jl"])
    @test "test_model.jl" in model && "test_physics.jl" in model && "test_golden.jl" ∉ model
    # Shared code, the module file and unknown paths run everything
    all_files = select_tests()
    @test tests_for(["test/support/testutils.jl"])[1] == all_files
    @test tests_for(["src/GREBClimate.jl"])[1] == all_files
    @test startswith(tests_for(["src/newdir/new.jl"])[2], "full: no rule")
    # A runner script has no test of its own beyond this file
    @test tests_for(["test/support/parallel.jl"]) == (["test_selection.jl"], "mapped")
    # Several paths: the union, in SHARD order
    @test tests_for(["src/io.jl", "src/core/calendar.jl"])[1] == ["test_calendar.jl", "test_io.jl"]
end
