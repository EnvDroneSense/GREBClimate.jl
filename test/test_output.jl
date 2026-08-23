# Per-timestep bookkeeping: diagnostics!, output!, time_loop!, climatology helpers.

@testset "build_monthly_climatology/apply_scenario_anomalies" begin
    # Hand-built records with every field filled to one scalar value
    # make the averaging arithmetic trivial to check by hand.
    mkrec(v) =(Ts=fill(v, GREBClimate.xdim, GREBClimate.ydim), Ta=fill(v, GREBClimate.xdim, GREBClimate.ydim),
        To=fill(v, GREBClimate.xdim, GREBClimate.ydim), q=fill(v, GREBClimate.xdim, GREBClimate.ydim),
        albedo=fill(v, GREBClimate.xdim, GREBClimate.ydim), ice=fill(v, GREBClimate.xdim, GREBClimate.ydim),
        precip=fill(v, GREBClimate.xdim, GREBClimate.ydim), evap=fill(v, GREBClimate.xdim, GREBClimate.ydim),
        qcrcl=fill(v, GREBClimate.xdim, GREBClimate.ydim), sw=fill(v, GREBClimate.xdim, GREBClimate.ydim),
        lw=fill(v, GREBClimate.xdim, GREBClimate.ydim), qlat=fill(v, GREBClimate.xdim, GREBClimate.ydim),
        qsens=fill(v, GREBClimate.xdim, GREBClimate.ydim))

    @test build_monthly_climatology(MonthlyRecord[]) == MonthlyRecord[]

    # Final-year-only.
    two_years = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
    clim = build_monthly_climatology(two_years)
    @test length(clim) == 12
    for m in 1:12
        @test all(==(Float64(m + 12)), clim[m].Ts)
    end

    # Non-12-multiple count: only 5 records -> months 6..12 never occur
    # and must fall back to records[1] exactly (not e.g. NaN/zero).
    five = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:5]
    clim5 = build_monthly_climatology(five)
    for m in 1:5
        @test all(==(Float64(m)), clim5[m].Ts)
    end
    for m in 6:12
        @test clim5[m] == five[1]
    end

    ctrl_clim = MonthlyRecord[mkrec(100.0 + m) for m in 1:12]
    scnr = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
    anom = apply_scenario_anomalies(scnr, ctrl_clim)
    @test all(==(1.0 - 101.0), anom[1].Ts)
    @test all(==(13.0 - 101.0), anom[13].Ts)
    @test all(==(12.0 - 112.0), anom[12].Ts)

    # Early-return guards: empty scnr_records or empty ctrl_clim ->
    # scnr_records passed straight through, not turned into anomalies.
    @test apply_scenario_anomalies(MonthlyRecord[], ctrl_clim) == MonthlyRecord[]
    @test apply_scenario_anomalies(scnr, MonthlyRecord[]) == scnr
end

@testset "compute_annual_ice_climatology" begin
    # Same record-index-encodes-value trick as the climatology test
    # above: final-year-only.
    mkrec(v) = (Ts=zeros(GREBClimate.xdim, GREBClimate.ydim), Ta=zeros(GREBClimate.xdim, GREBClimate.ydim),
        To=zeros(GREBClimate.xdim, GREBClimate.ydim), q=zeros(GREBClimate.xdim, GREBClimate.ydim),
        albedo=zeros(GREBClimate.xdim, GREBClimate.ydim), ice=fill(v, GREBClimate.xdim, GREBClimate.ydim),
        precip=zeros(GREBClimate.xdim, GREBClimate.ydim), evap=zeros(GREBClimate.xdim, GREBClimate.ydim),
        qcrcl=zeros(GREBClimate.xdim, GREBClimate.ydim), sw=zeros(GREBClimate.xdim, GREBClimate.ydim),
        lw=zeros(GREBClimate.xdim, GREBClimate.ydim), qlat=zeros(GREBClimate.xdim, GREBClimate.ydim),
        qsens=zeros(GREBClimate.xdim, GREBClimate.ydim))

    @test all(iszero, compute_annual_ice_climatology(MonthlyRecord[]))

    two_years = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
    clim = compute_annual_ice_climatology(two_years)
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
    diagnostics!(1, 1970, 340.0, surf, tend, fields, state, ts)
    @test all(==(280.0), state.Tsmn)   # accumulated once, no averaging/reset yet

    ts.ityr = GREBClimate.nstep_yr
    captured = mktemp() do path, io
        redirect_stdout(io) do
            diagnostics!(GREBClimate.nstep_yr, 1970, 340.0, surf, tend, fields, state, ts)
        end
        flush(io)
        read(path, String)
    end
    @test occursin("1970", captured)   # prints the annual summary line
    @test all(iszero, state.Tsmn)      # reset after year end
end

@testset "output! pushes a monthly-mean MonthlyRecord at month boundaries" begin
    ws = CirculationWorkspace()
    acc = MonthlyAccumulator()
    ts = TimeState(1, 1)
    surf = SurfaceState(fill(280.0, GREBClimate.xdim, GREBClimate.ydim), fill(270.0, GREBClimate.xdim, GREBClimate.ydim),
        fill(285.0, GREBClimate.xdim, GREBClimate.ydim), fill(0.005, GREBClimate.xdim, GREBClimate.ydim))
    tend = (albedo=fill(0.3, GREBClimate.xdim, GREBClimate.ydim), SW=fill(100.0, GREBClimate.xdim, GREBClimate.ydim),
        ice_cover=fill(0.1, GREBClimate.xdim, GREBClimate.ydim), LW_surf=fill(-50.0, GREBClimate.xdim, GREBClimate.ydim),
        Q_lat=fill(-20.0, GREBClimate.xdim, GREBClimate.ydim), Q_sens=fill(-5.0, GREBClimate.xdim, GREBClimate.ydim))
    ws.precip_out .= 2.0
    ws.evap_out .= 1.0
    ws.qcrcl_out .= 0.5

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
end

@testset "time_loop! integrates one timestep and clamps at min_T_K" begin
    fields = ClimateFields()
    fields.z_topo .= -1.0
    fields.mldclim .= 50.0
    fields.Tclim .= 280.0
    fields.Toclim .= 285.0
    fields.qclim .= 0.006
    fields.cldclim .= 0.5
    fields.swetclim .= 0.5
    fields.uclim .= 2.0
    fields.vclim .= 1.0
    fields.omegaclim .= 0.001
    fields.omegastdclim .= 0.01
    fields.wsclim .= 4.0
    cfg = create_experiment_config(:full_model)
    ini = init_model!(cfg, fields)

    state = ModelState()
    ws = CirculationWorkspace()
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
