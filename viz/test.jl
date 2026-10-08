# Tests for GREBViz on synthetic records; no dataset needed.
#
#   julia --project=viz viz/test.jl

using Test
include(joinpath(@__DIR__, "GREBViz.jl"))
using .GREBViz

const NX, NY = 96, 48
rec(f) = (Ts=Float32[f(i, j) for i in 1:NX, j in 1:NY], q=fill(1f-3, NX, NY))
fields = (z_topo=Float32[i <= NX ÷ 2 ? 100 : -100 for i in 1:NX, j in 1:NY],)

@testset "GREBViz" begin
    @testset "grid" begin
        @test lats(NY)[1] ≈ -90 + 180 / NY / 2
        @test lats(NY) ≈ -reverse(lats(NY))
        @test lons(NX)[end] ≈ 360 - 360 / NX / 2
        @test size(area_weights(NX, NY)) == (NX, NY) && all(>(0), area_weights(NX, NY))
    end

    @testset "reductions" begin
        @test all(≈(288), series([rec((i, j) -> 288) for _ in 1:5], :Ts))
        hot = [rec((i, j) -> 300 - abs(lats(NY)[j]))]                     # equator-hot field
        @test series(hot, :Ts)[1] > sum(hot[1].Ts) / length(hot[1].Ts)   # weighted > plain mean
        @test annual(collect(1.0:30)) == [6.5, 18.5]                      # partial year dropped
        @test seasonal_cycle([rec((i, j) -> m % 12) for m in 0:35], :Ts) ≈ 0:11
        @test isempty(seasonal_cycle([rec((i, j) -> 1) for _ in 1:11], :Ts))

        # With times, the month and the year come from them, not from the position:
        # 24 records from July 1950 hold one whole year, 1951 (records 7 to 18)
        times = [(year=1950 + (m + 5) ÷ 12, month=mod1(m + 6, 12)) for m in 1:24]
        @test GREBViz.by_position(13)[[1, 12, 13]] == [(year=1, month=1), (year=1, month=12), (year=2, month=1)]
        @test GREBViz.yearblocks(times) == [collect(7:18)]
        @test annual(collect(1.0:24), times) == [12.5]
        by_month = [rec((i, j) -> t.month) for t in times]
        @test seasonal_cycle(by_month, :Ts, times) ≈ 1:12
        @test length(map_frames(by_month, :Ts; step=:year, times=times)) == 1

        recs = [rec((i, j) -> t) for t in 1:4]
        @test all(==(2.5f0), field(recs, :Ts))
        @test all(==(4f0), field(recs, :Ts; month=:last))
        @test all(==(2f0), field(recs, :Ts; month=2))
        @test_throws BoundsError field(recs, :Ts; month=9)
        months, φ, m = hovmoller(recs, :Ts)
        @test size(m) == (4, NY) && m[:, 1] ≈ 1:4 && φ == lats(NY)
        @test length(map_frames([rec((i, j) -> 1) for _ in 1:30], :Ts; step=:year)) == 2
    end

    ctrl = [rec((i, j) -> 280 + m % 12) for m in 0:23]
    anom = [rec((i, j) -> 1f0 + m / 10) for m in 0:23]
    res = (ctrl=ctrl, scnr=anom)

    @testset "runs" begin
        @test [r[1] for r in GREBViz.runs(res)] == ["control", "scenario anomaly"]
        @test GREBViz.runs((ctrl=ctrl, scnr=ctrl))[2][3] == false          # orbital: stays absolute
        @test length(GREBViz.runs((ctrl=ctrl, scnr=typeof(ctrl)()))) == 1
        @test_throws ErrorException GREBViz.runs((ctrl=typeof(ctrl)(), scnr=typeof(ctrl)()))
        # A result with record times is dated; one without counts from year 1
        scnr_time = [(year=1950 + m ÷ 12, month=m % 12 + 1) for m in 0:23]
        dated = GREBViz.runs((ctrl=ctrl, scnr=anom, ctrl_time=GREBViz.by_position(24), scnr_time=scnr_time))
        @test dated[2][4] == scnr_time && dated[2][5]
        @test GREBViz.runs(res)[2][4] == GREBViz.by_position(24) && !GREBViz.runs(res)[2][5]
        @test evolution((ctrl=ctrl, scnr=anom, ctrl_time=GREBViz.by_position(24), scnr_time=scnr_time)).panels[2].labels[[1, 24]] ==
              ["Jan 1950", "Dec 1951"]
        @test evolution(res; step=:year).panels[1].labels == ["year 1", "year 2"]
        @test evolution(res).panels[1].labels[13] == "Jan, year 2"
    end

    @testset "plots" begin
        for x in (res, ctrl)
            @test plot_map(x; fields=fields) isa GREBViz.Plots.Plot
            @test plot_map(x; month=:last) isa GREBViz.Plots.Plot
            @test plot_timeseries(x; annual=true) isa GREBViz.Plots.Plot
            @test plot_seasonal(x; var=:q) isa GREBViz.Plots.Plot
            @test plot_hovmoller(x) isa GREBViz.Plots.Plot
        end
        ev = evolution(res; step=:year)
        @test frame_count(ev) == 2
        @test ev.panels[2].clims[1] == -ev.panels[2].clims[2]              # anomaly: symmetric
        @test evolution_frame(ev, 3; fields=fields) isa GREBViz.Plots.Plot
        @test_throws ErrorException evolution(ctrl[1:11]; step=:year)      # no whole year
    end
end
