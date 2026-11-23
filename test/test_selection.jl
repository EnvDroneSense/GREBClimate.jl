# The source-to-test map (support/selection.jl) that test/support/changed.jl reads.

@testset "every source file has a rule; every rule names real test files" begin
    root = normpath(joinpath(@__DIR__, ".."))
    sources = [replace(relpath(joinpath(dir, f), root), "\\" => "/")
               for top in ("src", "ext") for (dir, _, fs) in walkdir(joinpath(root, top)) for f in fs if endswith(f, ".jl")]
    @test !isempty(sources)
    named(path) = any(startswith(path, p) for p in FULL_PATHS) || any(r -> startswith(path, r[1]), SOURCE_MAP)
    @test filter(!named, sources) == String[]
    for (_, names) in SOURCE_MAP, name in names

        name == "@standard" || @test resolve_test_file(name) in first.(SHARD)
    end
    # Every path prefix names something that exists
    for prefix in [first.(SOURCE_MAP); FULL_PATHS[1:5]; RUNNER_PATHS]
        @test ispath(joinpath(root, prefix))
    end
end
