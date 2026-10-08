# The experiment helpers in tools/experiments/helpers.jl.

const ET = let tools = Module()
    Base.include(tools, joinpath(@__DIR__, "..", "tools", "experiments", "helpers.jl"))
    tools.ExperimentTools
end

@testset "gregory recovers a straight line" begin
    dT = collect(0.5:0.25:3.0)
    g = ET.gregory(dT, 5.0 .- 2.0 .* dT)
    @test g.forcing ≈ 5.0 && g.feedback ≈ -2.0 && g.sensitivity ≈ 2.5 && g.r2 ≈ 1
    @test g.reached == dT[end] && g.remaining ≈ 5.0 - 2.0 * dT[end]
    @test_throws ArgumentError ET.gregory([1.0, 2.0], [1.0, 2.0])
end

@testset "map diagnostics" begin
    ones_map = fill(1.0f0, X, Y)
    @test ET.area_mean(ones_map) ≈ 1
    @test ET.area_mean(ones_map) ≈ global_mean(ones_map)
    @test ET.zonal_mean(ones_map) == fill(1.0, Y)
    flat = ET.polar_amplification(ones_map)
    @test flat.north ≈ 1 && flat.south ≈ 1

    pole = [abs(GREBClimate.lat_grid[j]) >= 60 ? 2.0f0 : 1.0f0 for _ in 1:X, j in 1:Y]
    amp = ET.polar_amplification(pole)
    @test amp.north ≈ amp.south > 1
    @test amp.north ≈ 2 / ET.area_mean(pole)

    z_topo = synthetic_fields().z_topo
    land = Float32.(GREBClimate.is_land.(z_topo))
    @test any(land .== 1) && any(land .== 0)
    lo = ET.land_ocean(land, z_topo)
    @test lo.land ≈ 1 && abs(lo.ocean) < 1e-12
    @test ET.area_mean(land, land .== 1) ≈ 1
    @test_throws ArgumentError ET.area_mean(land, falses(X, Y))

    a = [Float32(sin(i / 7) + cos(j / 5)) for i in 1:X, j in 1:Y]
    @test ET.pattern_correlation(a, a) ≈ 1
    @test ET.pattern_correlation(a, 3 .* a .+ 2) ≈ 1
    @test ET.pattern_correlation(a, -a) ≈ -1
    @test_throws DimensionMismatch ET.area_mean(zeros(Float32, 3, 3))
end

@testset "sine_fit finds amplitude, lag and mean" begin
    period, lag = 40, 5.0
    series = [2 + 3 * sin(2π * (t - lag) / period) for t in 0:159]
    fit = ET.sine_fit(series, period)
    @test fit.amplitude ≈ 3 && fit.lag ≈ lag && fit.mean ≈ 2
end

@testset "annual means group by record time" begin
    times = [(year = 1950 + (m - 1) ÷ 12, month = mod1(m, 12)) for m in 1:30]
    records = [uniform_record(Float32(m)) for m in 1:30]
    years, values = ET.annual_means(records, times, :Ts)
    @test years == [1950, 1951]                      # 1952 has six months only
    @test values ≈ [gmean(1:12), gmean(13:24)]
    @test ET.annual_means(records, times, rec -> 2 * rec.Ts[1, 1])[2] ≈ 2 .* values
    @test_throws DimensionMismatch ET.annual_means(records, times[1:5], :Ts)
end

@testset "in_range_report" begin
    times = [(year = 1950, month = m) for m in 1:3]
    fine = [uniform_record(280f0; q = 0.01f0) for _ in 1:3]
    @test ET.in_range_report(fine, times, false) === nothing
    runaway = [fine[1], uniform_record(280f0; q = 0.01f0, Ts = 1f24), fine[3]]
    @test ET.in_range_report(runaway, times, false) == (year = 1950, month = 2, field = :Ts, value = 1f24)
    nan = [fine[1], fine[2], uniform_record(280f0; q = 0.01f0, Ta = NaN32)]
    @test ET.in_range_report(nan, times, false).month == 3
    # An anomaly of a few kelvin is fine; the absolute range does not apply to it
    small = [uniform_record(2f0; q = 0.001f0) for _ in 1:3]
    @test ET.in_range_report(small, times, true) === nothing
    @test ET.in_range_report(small, times, false).field == :Ts
    @test ET.in_range_report([uniform_record(1f24; q = 0.001f0)], times[1:1], true).field == :Ts
end

@testset "summarize reduces a model result" begin
    result = run_synthetic(RunSpec(ctrl = 1, scnr = 2), preset(:co2_double; corrections = NoCorrections()))
    s = ET.summarize(result)
    @test s.scnr_anomaly && s.ctrl.years == [1970] && s.scnr.years == [1950, 1951]
    @test all(length(s.scnr[k]) == 2 for k in (:Ts, :precip, :ice, :net))
    @test s.ctrl.Ts[1] ≈ gmean([global_mean(r.Ts) for r in result.ctrl])
    @test s.scnr.net[2] ≈ gmean([global_mean(r.sw) - global_mean(r.olr) for r in result.scnr[13:24]])
    @test s.out_of_range == (ctrl = nothing, scnr = nothing)

    # A scenario can leave the range while its control is fine
    broken = merge(result, (scnr = [uniform_record(1f24; q = 0.001f0); result.scnr[2:end]],))
    @test ET.summarize(broken).out_of_range.ctrl === nothing
    @test ET.summarize(broken).out_of_range.scnr.field == :Ts
end
