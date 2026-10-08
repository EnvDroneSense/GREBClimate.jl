# Per-timestep bookkeeping: diagnostics!, output!, time_loop!, climatology helpers.

@testset "monthly_climatology/scenario_anomalies" begin
    # Records with every field filled to one scalar value make the averaging
    # arithmetic trivial to check by hand.
    mkrec = uniform_record

    @test monthly_climatology(MonthlyRecord[]) == MonthlyRecord[]

    # Final-year-only.
    two_years = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
    clim = monthly_climatology(two_years)
    @test length(clim) == 12
    for m in 1:12
        @test all(==(Float64(m + 12)), clim[m].Ts)
    end

    # Non-12-multiple count: only 5 records -> months 6..12 never occur
    # and must fall back to records[1] exactly (not e.g. NaN/zero).
    five = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:5]
    clim5 = monthly_climatology(five)
    for m in 1:5
        @test all(==(Float64(m)), clim5[m].Ts)
    end
    for m in 6:12
        @test clim5[m] == five[1]
    end

    ctrl_clim = MonthlyRecord[mkrec(100.0 + m) for m in 1:12]
    scnr = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
    anom = scenario_anomalies(scnr, ctrl_clim)
    @test all(==(1.0 - 101.0), anom[1].Ts)
    @test all(==(13.0 - 101.0), anom[13].Ts)
    @test all(==(12.0 - 112.0), anom[12].Ts)

    # Early-return guards: empty scnr_records or empty ctrl_clim ->
    # scnr_records passed straight through, not turned into anomalies.
    @test scenario_anomalies(MonthlyRecord[], ctrl_clim) == MonthlyRecord[]
    @test scenario_anomalies(scnr, MonthlyRecord[]) == scnr
end

@testset "ice_climatology" begin
    # Same record-index-encodes-value trick as the climatology test
    # above: final-year-only.
    mkrec(v) = uniform_record(0; ice = v)

    @test all(iszero, ice_climatology(MonthlyRecord[]))

    two_years = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
    clim = ice_climatology(two_years)
    @test size(clim) == (GREBClimate.xdim, GREBClimate.ydim, 12)
    for m in 1:12
        @test all(==(Float64(m + 12)), clim[:, :, m])
    end
end

@testset "diagnostics! accumulates annual means and resets at year end" begin
    fields = ClimateFields()
    state = ModelState()
    ts = TimeState(1, 1)
    z = () -> zeros(GREBClimate.xdim, GREBClimate.ydim)
    surf = SurfaceState(fill(280.0, GREBClimate.xdim, GREBClimate.ydim), fill(270.0, GREBClimate.xdim, GREBClimate.ydim),
        fill(285.0, GREBClimate.xdim, GREBClimate.ydim), fill(0.005, GREBClimate.xdim, GREBClimate.ydim))
    tend = (albedo=fill(0.3, GREBClimate.xdim, GREBClimate.ydim), SW=fill(100.0, GREBClimate.xdim, GREBClimate.ydim),
        ice_cover=z(), LW_surf=fill(-50.0, GREBClimate.xdim, GREBClimate.ydim),
        Q_lat=fill(-20.0, GREBClimate.xdim, GREBClimate.ydim), Q_sens=fill(-5.0, GREBClimate.xdim, GREBClimate.ydim),
        Q_lat_air=fill(20.0, GREBClimate.xdim, GREBClimate.ydim), dq_eva=z(),
        dq_rain=z(), dq_crcl=z(), dTa_crcl=z(), dT_ocean=z(), dTo=z(),
        LW_down=fill(30.0, GREBClimate.xdim, GREBClimate.ydim), LW_up=fill(80.0, GREBClimate.xdim, GREBClimate.ydim),
        em=fill(0.9, GREBClimate.xdim, GREBClimate.ydim))

    ts.ityr = 1
    diagnostics!(1970, surf, state, ts)
    @test all(==(280.0), state.Ts_annual_mean)   # accumulated once, no averaging/reset yet

    ts.ityr = GREBClimate.nstep_yr
    # Logs the annual summary line. Two steps of 280 K over a year of 730: 0.77 K,
    # the same in the global mean and in the two cells
    @test_logs (:info, "1970: Ts global mean -272.38 °C; 178 E 9 N -272.38; 58 E 51 N -272.38") diagnostics!(
        1970, surf, state, ts)
    @test all(iszero, state.Ts_annual_mean)      # reset after year end
end

@testset "output! pushes a monthly-mean MonthlyRecord at month boundaries" begin
    ws = ModelWorkspace()
    acc = MonthlyAccumulator()
    ts = TimeState(1, 1)
    surf = SurfaceState(fill(280.0, GREBClimate.xdim, GREBClimate.ydim), fill(270.0, GREBClimate.xdim, GREBClimate.ydim),
        fill(285.0, GREBClimate.xdim, GREBClimate.ydim), fill(0.005, GREBClimate.xdim, GREBClimate.ydim))
    tend = (albedo=fill(0.3, GREBClimate.xdim, GREBClimate.ydim), SW=fill(100.0, GREBClimate.xdim, GREBClimate.ydim),
        ice_cover=fill(0.1, GREBClimate.xdim, GREBClimate.ydim), LW_surf=fill(-50.0, GREBClimate.xdim, GREBClimate.ydim),
        Q_lat=fill(-20.0, GREBClimate.xdim, GREBClimate.ydim), Q_sens=fill(-5.0, GREBClimate.xdim, GREBClimate.ydim),
        LW_up=fill(-200.0, GREBClimate.xdim, GREBClimate.ydim), LW_down=fill(-210.0, GREBClimate.xdim, GREBClimate.ydim),
        em=fill(0.75, GREBClimate.xdim, GREBClimate.ydim))
    ws.precip .= 2.0
    ws.evap .= 1.0
    ws.qcrcl .= 0.5

    output_buf = MonthlyRecord[]
    irec, mon = 0, 1
    ndt = GREBClimate.ndt_days
    ndays_jan = GREBClimate.cjday_mon[1]
    for day in 1:ndays_jan, step in 1:ndt
        it = (day - 1) * ndt + step
        ts.jday = day
        (mon, irec) = output!(it, irec, mon, surf, tend, ws, output_buf, acc, ts)
    end

    @test length(output_buf) == 1
    @test irec == 1
    @test mon == 2
    @test all(==(280.0), output_buf[1].Ts)
    @test all(==(2.0), output_buf[1].precip)
    # Out to space: what the air emits upward plus the part of the surface's
    # emission the air lets through, 200 + (1 - 0.75) * 50. Positive upward.
    @test all(==(212.5), output_buf[1].olr)
    @test all(==(210.0), output_buf[1].lwdown)   # positive into the surface

    # Every field of the list reaches the record under its own name
    @test keys(output_buf[1]) == GREBClimate.OUTPUT_FIELDS
    expected = (Ts=280.0, Ta=270.0, To=285.0, q=0.005, albedo=0.3, ice=0.1, precip=2.0, evap=1.0,
        qcrcl=0.5, sw=100.0, lw=-50.0, qlat=-20.0, qsens=-5.0, olr=212.5, lwdown=210.0)
    for name in GREBClimate.OUTPUT_FIELDS
        @test all(x -> isapprox(x, expected[name]; rtol=1.0f-5), output_buf[1][name])
    end
    @test all(f -> all(iszero, getfield(acc, f)), fieldnames(MonthlyAccumulator))   # reset after the record
end

# All ocean, with wind and rising air
_time_loop_fields() = constant_fields(z_topo = -1.0, swet = 0.5, u = 2.0, v = 1.0, omega = 0.001, omega_std = 0.01, ws = 4.0)

@testset "time_loop! integrates one timestep and clamps at min_T_K" begin
    fields = _time_loop_fields()
    cfg = resolve(preset(:full_model))
    ini = init_model!(cfg, fields)

    state = ModelState()
    ws = ModelWorkspace()
    acc = MonthlyAccumulator()
    ts = TimeState(1, 1)

    Ts = fill(GREBClimate.min_T_K - 0.5, GREBClimate.xdim, GREBClimate.ydim)
    Ta = copy(ini.Ta_ini)
    To = copy(ini.To_ini)
    q = copy(ini.q_ini)
    output_buf = MonthlyRecord[]

    (mon, irec) = time_loop!(1, 1970, ini.CO2_ctrl, 1, 0, Ts, Ta, q, To, output_buf,
        fields, state, ws, acc, ts, cfg)

    @test all(isfinite, Ts)
    @test all(isfinite, Ta)
    @test all(isfinite, To)
    @test all(isfinite, q)
    @test all(>=(GREBClimate.min_T_K), Ts)
    @test all(>=(GREBClimate.min_T_K), Ta)
    @test mon == 1
    @test irec == 0
end

@testset "time_loop!'s min_T_K floor leaves NaN as NaN" begin
    # A max-based floor inside @turbo turns NaN into min_T_K, which hides a
    # failed run behind a plausible-looking cold planet.
    fields = _time_loop_fields()
    cfg = resolve(preset(:full_model))
    ini = init_model!(cfg, fields)
    Ts = copy(ini.Ts_ini)
    Ta = copy(ini.Ta_ini)
    Ts[5, 5] = NaN32
    Ta[40, 30] = NaN32

    time_loop!(1, 1970, ini.CO2_ctrl, 1, 0, Ts, Ta, copy(ini.q_ini), copy(ini.To_ini),
        MonthlyRecord[], fields, ModelState(), ModelWorkspace(), MonthlyAccumulator(),
        TimeState(1, 1), cfg)

    @test isnan(Ts[5, 5])
    @test isnan(Ta[40, 30])
end

@testset "global_mean weights each row by the cosine of its latitude" begin
    X, Y = GREBClimate.xdim, GREBClimate.ydim
    @test global_mean(fill(3.5f0, X, Y)) ≈ 3.5
    # A field equal to the latitude weight: mean of cos^2 over mean of cos
    w = cosd.(GREBClimate.lat_grid)
    field = repeat(w', X, 1)
    @test global_mean(field) ≈ sum(abs2, w) / sum(w) rtol = 1e-6
    # The poles count less than in a plain mean
    polar = zeros(Float32, X, Y); polar[:, 1] .= 1; polar[:, Y] .= 1
    @test global_mean(polar) < sum(polar) / length(polar) / 10
    @test_throws DimensionMismatch global_mean(zeros(Y, X))
    # The sample cells of the annual line are where their labels say
    for (label, i, j) in GREBClimate._SAMPLE_CELLS
        lon, lat = (i - 0.5) * GREBClimate.dlon, GREBClimate.lat_grid[j]
        @test label == "$(round(Int, lon)) E $(round(Int, lat)) N"
    end
end
