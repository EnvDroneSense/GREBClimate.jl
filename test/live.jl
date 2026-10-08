# Iterate on tests in one long-lived Julia session, so the package and the test
# files compile once instead of on every run:
#
#   julia --project=. -i test/live.jl
#   julia> t("state")            # one file; name as in SHARD, prefix and .jl optional
#   julia> t("state", "io")      # several
#   julia> t(:smoke)             # a tier: :smoke, :standard or :full
#
# Revise (if installed in your global environment) picks up edits to src/. Restart
# after changing a struct definition. This session does not force
# --check-bounds=yes the way Pkg.test() does, so an out-of-bounds index can pass
# here and fail there: run the full Pkg.test() before committing.

try
    @eval using Revise
catch
    @info "Revise not available; restart Julia after editing src/"
end

using GREBClimate
using GREBClimate: SWradiation!, LWradiation!, hydro!, convergence!, seaice!, deep_ocean!,
    diffusion!, advection!, circulation!, tendencies!, time_loop!, output!, diagnostics!,
    qflux_correction!, init_model!
using Test

include(joinpath(@__DIR__, "support", "testutils.jl"))
include(joinpath(@__DIR__, "support", "selection.jl"))

"""Run the named test files, or a whole tier (`t(:smoke)`), in this session."""
function t(names::AbstractString...)
    files = select_tests(; names = collect(names))
    return _run(files)
end
t(tier::Symbol) = _run(select_tests(; tier = String(tier)))

function _run(files)
    isempty(files) && return nothing
    elapsed = @elapsed @testset "live" begin
        for file in files
            @testset "$file" begin
                include(joinpath(@__DIR__, file))
            end
        end
    end
    println("ran ", join(files, ", "), " in ", round(elapsed; digits = 1), " s (subset only, not the gate)")
    return nothing
end
