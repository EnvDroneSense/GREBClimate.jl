# Enzyme rules for fastlog and fastexp (ext/GREBEnzymeExt.jl). Enzyme is not a
# test dependency: the check runs in a subprocess when Enzyme can be found, in
# the default environment or in a shared environment named `enzyme`
# (`julia --project=@enzyme -e 'using Pkg; Pkg.add("Enzyme")'`), and is skipped
# otherwise. Neither CI nor a fresh checkout has Enzyme, so there it is skipped.

enzyme_env = joinpath(first(DEPOT_PATH), "environments", "enzyme")

@testset "fastlog/fastexp: Enzyme rules match the analytic derivatives" begin
    if Base.find_package("Enzyme") === nothing && !isdir(enzyme_env)
        @info "Enzyme not installed; skipping test_enzyme.jl"
        @test_skip true
    else
        script = """
        using GREBClimate, Enzyme
        using GREBClimate: fastlog, fastexp
        sum_log(x) = (s = 0.0f0; for i in eachindex(x); s += fastlog(x[i]); end; s)
        sum_exp(x) = (s = 0.0f0; for i in eachindex(x); s += fastexp(x[i]); end; s)
        x = Float32[0.5, 2.0, 7.0]
        println(Enzyme.autodiff(Forward, fastlog, Duplicated, Duplicated(2.0f0, 1.0f0))[1] ≈ 0.5f0)
        println(Enzyme.autodiff(Forward, fastexp, Duplicated, Duplicated(0.3f0, 1.0f0))[1] ≈ exp(0.3f0))
        println(Enzyme.autodiff(Reverse, fastlog, Active, Active(2.0f0))[1][1] ≈ 0.5f0)
        println(Enzyme.autodiff(Reverse, fastexp, Active, Active(0.3f0))[1][1] ≈ exp(0.3f0))
        for (f, df) in ((sum_log, inv), (sum_exp, exp))
            dx = zero(x)
            Enzyme.autodiff(Reverse, f, Active, Duplicated(copy(x), dx))
            println(dx ≈ df.(x))
        end
        """
        exe = first(Base.julia_cmd())
        project = normpath(joinpath(@__DIR__, ".."))
        load_path = join(["@", enzyme_env, "@v#.#", "@stdlib"], Sys.iswindows() ? ";" : ":")
        cmd = addenv(`$exe --startup-file=no --project=$project -e $script`, "JULIA_LOAD_PATH" => load_path)
        out = read(cmd, String)
        @test split(out) == fill("true", 6)
    end
end
