# Every experiment preset still imposes what the reference recorded
# (test/data/preset_reference.jl, written by tools/validation/preset_reference.jl).

include("presetcapture.jl")
include(joinpath("data", "preset_reference.jl"))

@testset "preset reference: $p" for p in PRESETS
    got, ref = capture_preset(p), PRESET_REFERENCE[p]
    for k in (:co2_ctrl, :start_year, :output, :co2_table, :solar_table, :boundary, :static_mask, :dynamic_mask)
        @test getfield(got, k) == getfield(ref, k)
    end
    # Float32 sin/cos may differ in the last digit across Julia versions
    @test isapprox(unrle(got.co2), unrle(ref.co2); rtol = 1e-6)
    @test isapprox(unrle(got.solar), unrle(ref.solar); rtol = 1e-6)
end

@testset "preset reference covers every preset" begin
    @test Set(PRESETS) == Set(keys(PRESET_REFERENCE))
end
