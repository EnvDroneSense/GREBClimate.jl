# A run on a synthetic fixture against stored reference values. Runs in CI, which has no
# dataset, so a change to a kernel term that only the golden test sees still fails here.

include(joinpath(@__DIR__, "support", "regression_reference.jl"))

@testset "spin-up, control and 2xCO2 run on the rich fixture keep their reference means" begin
    res = run_synthetic(RunSpec(ctrl=1, scnr=1), preset(:co2_double; corrections=SpinUp(1)); fields=rich_fields())
    records = vcat(res.ctrl, res.scnr)
    @test length(records) == 24 && res.scnr_anomaly
    # The bound is a share of the field's largest reference value. The check-bounds and
    # plain builds differ by 2e-4 of it at most, except for `qcrcl`, a small difference of
    # large terms that moves by 2 %; a change in a kernel term moves a field by more.
    for name in fieldnames(MonthlyRecord)
        ref = getfield(REGRESSION_REFERENCE, name)
        got = [gmean(getfield(records[i], name)) for i in REGRESSION_RECORDS]
        share = name === :qcrcl ? 0.05 : 1e-3
        @test all(abs.(got .- ref) .<= share * maximum(abs, ref))
    end
end
