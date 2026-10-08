using GREBClimate
# The kernels and loop functions are internal: not exported, reached by name
using GREBClimate: SWradiation!, LWradiation!, hydro!, convergence!, seaice!, deep_ocean!,
    diffusion!, advection!, circulation!, tendencies!, time_loop!, output!, diagnostics!,
    qflux_correction!, init_model!
using Test

include(joinpath("support", "testutils.jl"))
include(joinpath("support", "selection.jl"))

# One file per subject; `SHARD` (support/selection.jl) lists them. A file not
# listed there is never run. Select a subset, most specific first:
#   Pkg.test(test_args = ["state", "io"])  or  GREB_TEST_FILES=state,io
#   GREB_TEST_TIER=smoke|standard|full     (cumulative)
#   GREB_TEST_SHARD=light|heavy            (what CI runs)
# Unset means everything, which is the pre-commit gate. A subset run proves
# only the files it ran.
selected = select_tests(;
    names = isempty(ARGS) ? filter(!isempty, split(get(ENV, "GREB_TEST_FILES", ""), ',')) : ARGS,
    tier = get(ENV, "GREB_TEST_TIER", ""),
    shard = get(ENV, "GREB_TEST_SHARD", "all"),
)

@testset "GREBClimate.jl" begin
    for file in selected
        include(file)
    end
end
