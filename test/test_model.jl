# greb_model! integration: experiment dispatch, scenario tables, flux correction.

@testset "greb_model! baseline: default config runs to completion with the right output shape" begin
    cfg = create_experiment_config(:full_model)
    result = quiet() do
        greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = "", allow_uninitialized = true)
    end
    @test length(result.ctrl) == 12
    @test length(result.scnr) == 0
    @test result.ctrl[1] isa MonthlyRecord
end

@testset "greb_model! flux-correction spin-up: loaded files aren't overwritten by qflux_correction!" begin
    tmpdir = mktempdir()
    try
        mkpath(joinpath(tmpdir, "climatology"))
        GREBClimate.jldopen(joinpath(tmpdir, "climatology", "flux_corrections.jld2"), "w") do f
            f["Tsurf_flux_correction"] = fill(42.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
            f["vapour_flux_correction"] = fill(43.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
            f["Tocean_flux_correction"] = fill(44.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
        end

        cfg = create_experiment_config(:full_model)
        cfg.log_topo_drsp = false
        cfg.log_qflux_dmc = true
        fields = ClimateFields()
        quiet() do
            greb_model!(RunSpec(flux = 1, ctrl = 0, scnr = 0), cfg; jld2_dir = tmpdir, fields = fields, allow_uninitialized = true)
        end

        @test all(==(42.0), fields.TF_correct)
        @test all(==(43.0), fields.qF_correct)
        @test all(==(44.0), fields.ToF_correct)
    finally
        rm(tmpdir; recursive = true, force = true)
    end
end

@testset "greb_model! runs across log_eva / log_rain branches" begin
    # One end-to-end run proves the plumbing. The remaining log_eva values are
    # separate `@turbo` blocks in `hydro!`, reachable directly for the price of
    # one call each rather than a simulated year each.
    cfg = create_experiment_config(:full_model)
    cfg.log_eva, cfg.log_rain = -1, 0
    result = quiet() do
        greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = "", allow_uninitialized = true)
    end
    @test length(result.ctrl) == 12

    fields = synthetic_fields()
    Ts = fill(288.0f0, X, Y)
    q = fill(0.005f0, X, Y)
    for (log_eva, log_rain) in ((-1, -1), (0, 1), (1, 2), (2, 3))
        c = create_experiment_config(:full_model)
        c.log_eva, c.log_rain = log_eva, log_rain
        quiet() do
            init_model!(c, fields)
        end
        out = hydro!(Ts, q, fields, TimeState(1, 1), c, CirculationWorkspace())
        @test all(isfinite, out.Q_lat)
        @test all(isfinite, out.dq_rain)
    end
end

@testset "qflux_correction! pulls Ts/To/q to climatology; Ta gets no correction (matches Fortran)" begin
    fields = ClimateFields()
    fields.cap_surf .= GREBClimate.cap_ocean
    for j in 1:GREBClimate.ydim, i in 1:GREBClimate.xdim
        fields.Tclim[i, j, :] .= 280.0 + 5.0 * sin(i / 10.0) * cos(j / 8.0)
        fields.Toclim[i, j, :] .= 279.0
        fields.qclim[i, j, :] .= 0.006
    end
    cfg = create_experiment_config(:full_model)
    state = ModelState()
    ts = TimeState(1, 1)
    ws = CirculationWorkspace()

    Ts = fill(290.0, GREBClimate.xdim, GREBClimate.ydim)
    Ta = fill(290.0, GREBClimate.xdim, GREBClimate.ydim)
    q = fill(0.010, GREBClimate.xdim, GREBClimate.ydim)
    To = fill(285.0, GREBClimate.xdim, GREBClimate.ydim)

    GREBClimate.qflux_correction!(340.0, Ts, Ta, q, To, fields, state, ts, cfg, ws, 1)

    @test any(!=(0.0), fields.TF_correct)
    @test any(!=(0.0), fields.ToF_correct)
    @test any(!=(0.0), fields.qF_correct)
    @test all(isfinite, fields.TF_correct)
    @test all(isfinite, fields.ToF_correct)
    @test all(isfinite, fields.qF_correct)

    @test Ts ≈ fields.Tclim[:, :, 1]
    @test To ≈ fields.Toclim[:, :, 1]
    @test q ≈ fields.qclim[:, :, 1]

    @test all(isfinite, Ta)
    @test Ta != fill(290.0, GREBClimate.xdim, GREBClimate.ydim)
end

@testset "greb_model! swaps sw_solar for paleo experiments, restores after" begin
    # A paleo run's swapped solar table must not leak into a later run.
    fields = ClimateFields()
    saved_sw_solar = copy(fields.sw_solar)
    tmpdir = mktempdir()
    try
        mkpath(joinpath(tmpdir, "solar_scenarios"))
        distinctive_value = 999.0
        GREBClimate.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_paleo.jld2"), "w") do file
            file["data"] = fill(distinctive_value, GREBClimate.ydim, GREBClimate.nstep_yr)
            file["dim_names"] = ["lat", "time"]
        end

        cfg = create_experiment_config(:paleo_231kyr)
        captured = mktemp() do path, io
            result = redirect_stdout(io) do
                greb_model!(RunSpec(ctrl = 0), cfg; jld2_dir = tmpdir, fields = fields, allow_uninitialized = true)
            end
            flush(io)
            (result = result, text = read(path, String))
        end
        @test length(captured.result.scnr) == 12
        @test occursin("loading alternate solar forcing", captured.text)

        # sw_solar restored to its pre-run value after greb_model! returns
        @test fields.sw_solar == saved_sw_solar

        cfg_plain = create_experiment_config(:full_model)
        quiet() do
            greb_model!(RunSpec(scnr = 0), cfg_plain; jld2_dir = "", fields = fields, allow_uninitialized = true)
        end
        @test fields.sw_solar == saved_sw_solar
    finally
        rm(tmpdir; recursive = true, force = true)
    end
end

@testset "apply_dynamic_co2_mask! uses the annual-mean ice cover, not January" begin
    fields = ClimateFields()  # z_topo defaults to 0 everywhere -> land branch never fires
    cfg = PhysicsConfig(experiment = :regional_co2_ocean)

    icmn_ctrl = zeros(Float64, GREBClimate.xdim, GREBClimate.ydim, 12)
    # Cell A: January alone >= 0.5, but the other 11 months are 0 ->
    # annual mean ~0.083, NOT ice under the Fortran-matching rule.
    icmn_ctrl[1, 1, 1] = 1.0
    # Cell B: January alone < 0.5, but the other 11 months are 1.0 ->
    # annual mean ~0.917, IS ice under the Fortran-matching rule.
    icmn_ctrl[2, 1, 1] = 0.0
    icmn_ctrl[2, 1, 2:12] .= 1.0

    GREBClimate.apply_dynamic_co2_mask!(cfg, fields, icmn_ctrl)

    @test fields.co2_part[1, 1] == 1.0  # January said "ice"; annual mean says no
    @test fields.co2_part[2, 1] == 0.5  # January said "no ice"; annual mean says yes

    # forcing() no longer mutates the mask - it is pure for these two now.
    fresh = ClimateFields()
    forcing(1, 1970, cfg, fresh, icmn_ctrl)
    @test all(isone, fresh.co2_part)
    @test forcing(1, 1970, cfg, fresh, icmn_ctrl).CO2 == 680.0f0

    # The land/ice variant inverts the ocean mask and exempts ice cells.
    f_li = ClimateFields()
    cfg_li = PhysicsConfig(experiment = :regional_co2_land_ice)
    GREBClimate.apply_dynamic_co2_mask!(cfg_li, f_li, icmn_ctrl)
    @test f_li.co2_part[1, 1] == 0.5  # ocean cell, not annual-mean ice
    @test f_li.co2_part[2, 1] == 1.0  # annual-mean ice -> exempted back to 1.0

    # Every other experiment is a no-op.
    f_noop = ClimateFields()
    GREBClimate.apply_dynamic_co2_mask!(PhysicsConfig(experiment = :full_model),
                                        f_noop, icmn_ctrl)
    @test all(isone, f_noop.co2_part)
end

@testset "forcing()/init_model! dispatch on every experiment symbol" begin
    # These reach forcing() directly - no model run needed to prove the branch
    # exists and returns finite numbers.
    direct_dispatch_symbols = (
        :a1b_scenario, :co2_10x, :co2_half, :co2_zero, :solar_cycle_11yr,
        :co2_sine_wave, :co2_step, :modern_solar_paleo_co2,
        :earth_sun_distance, :regional_co2_nh, :regional_co2_sh,
        :regional_co2_tropics, :regional_co2_extratropics,
        :regional_co2_ocean, :regional_co2_land_ice, :regional_co2_winter,
        :regional_co2_summer,
    )
    for sym in direct_dispatch_symbols
        fields = ClimateFields()
        cfg = PhysicsConfig(experiment = sym)
        icmn_ctrl = zeros(Float64, X, Y, 1)
        quiet() do
            init_model!(cfg, fields)
        end
        result = forcing(1, 1970, cfg, fields, icmn_ctrl)
        @test isfinite(result.CO2)
        @test isfinite(result.sw_solar_forcing)
    end
end

@testset "IPCC scenario CO2 tables load under the right on-disk key" begin
    # `greb_model!` calls load_co2_scenario_jld2 once at scenario start. Assert
    # the loader against every key instead of paying a simulated year each;
    # :rcp60's on-disk key is "rcp6", and :historical_co2's is "hist".
    expected = Dict(:rcp26 => 400.0, :rcp45 => 401.0, :rcp60 => 402.0,
                    :ssp119 => 300.0, :ssp126 => 301.0, :ssp245 => 302.0,
                    :ssp460 => 303.0, :ssp585 => 304.0, :historical_co2 => 280.73)
    disk_key = Dict(:rcp60 => "rcp6", :historical_co2 => "hist")
    with_tempdir() do dir
        write_ipcc_scenarios(dir, Dict(
            get(disk_key, sym, string(sym)) => Dict(1950 => co2)
            for (sym, co2) in expected))
        for (sym, co2) in expected
            key = GREBClimate._CO2_SCENARIO_KEY
            scenario_key = get(key, sym, sym)
            table = load_co2_scenario_jld2(dir, scenario_key)
            @test isapprox(table[1950], co2; atol = 1e-3)
        end

        # One end-to-end run proves greb_model! actually wires the table in.
        cfg = PhysicsConfig(experiment = :rcp45)
        result = quiet() do
            greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = dir,
                        allow_uninitialized = true)
        end
        @test length(result.scnr) == 12
        @test cfg.co2_scenario == Dict(1950 => expected[:rcp45])

        # A year missing from the table must raise a clear error rather than
        # silently defaulting.
        cfg_missing_year = PhysicsConfig(experiment = :ssp585)
        @test_throws ErrorException greb_model!(RunSpec(ctrl = 0, scnr = 2),
            cfg_missing_year; jld2_dir = dir, allow_uninitialized = true)
    end
end

@testset "custom CO2 trajectory loads from a plain-text file" begin
    with_tempdir() do dir
        co2_path = joinpath(dir, "my_co2.txt")
        write(co2_path, "# comment line, should be skipped\n1950 300.0\n1951 301.0\n\n")

        # The parser is the thing under test - assert it directly.
        @test load_custom_co2_scenario(co2_path) == Dict(1950 => 300.0, 1951 => 301.0)

        # ...then one run to prove greb_model! reads cfg.custom_co2_path.
        cfg = create_experiment_config(:custom_co2; co2_path = co2_path)
        result = quiet() do
            greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = "",
                        allow_uninitialized = true)
        end
        @test length(result.scnr) == 12
        @test cfg.co2_scenario == Dict(1950 => 300.0, 1951 => 301.0)

        # Unset custom_co2_path must raise a clear error, not silently dispatch
        # or default.
        cfg_unset = PhysicsConfig(experiment = :custom_co2)
        @test_throws ErrorException greb_model!(RunSpec(ctrl = 0, scnr = 1),
            cfg_unset; jld2_dir = "", allow_uninitialized = true)
    end
end

@testset "paleo/orbital solar tables load; one runs end to end" begin
    with_tempdir() do dir
        write_solar_scenarios(dir)
        # The loader is the mechanism; assert all three tables directly.
        for (sym, ftype) in ((:paleo_solar_modern_co2, :paleo),
                             (:obliquity, :obliquity), (:eccentricity, :eccentricity))
            table = load_solar_forcing_jld2(dir, ftype, 0)
            @test size(table) == (Y, N)
            @test all(==(999.0f0), table)
        end

        cfg = PhysicsConfig(experiment = :obliquity)
        result = quiet() do
            greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = dir,
                        allow_uninitialized = true)
        end
        @test length(result.scnr) == 12
    end
end

@testset "greb_model! smoke: sst_plus1, decon presets, dynamic regional mask" begin
    cfg = PhysicsConfig(experiment = :sst_plus1)
    result = quiet() do
        greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = "",
                    allow_uninitialized = true)
    end
    @test length(result.scnr) == 12

    # :decon_mean_climate (control-run only) / :decon_2xco2 (scenario run).
    result_dmc = quiet() do
        greb_model!(RunSpec(ctrl = 1, scnr = 0),
                    create_experiment_config(:decon_mean_climate);
                    jld2_dir = "", allow_uninitialized = true)
    end
    @test length(result_dmc.ctrl) == 12

    result_drsp = quiet() do
        greb_model!(RunSpec(ctrl = 0, scnr = 1),
                    create_experiment_config(:decon_2xco2);
                    jld2_dir = "", allow_uninitialized = true)
    end
    @test length(result_drsp.scnr) == 12

    run_mask(sym) = begin
        f = ClimateFields()
        f.z_topo[1:(X - 48), :] .= 100.0f0   # left half land, right half ocean
        quiet() do
            greb_model!(RunSpec(ctrl = 1, scnr = 0), PhysicsConfig(experiment = sym);
                        jld2_dir = "", fields = f, allow_uninitialized = true)
        end
        f.co2_part
    end
    @test all(==(0.5f0), run_mask(:regional_co2_ocean))
    @test all(isone, run_mask(:regional_co2_land_ice))
end

@testset "CMIP5/ERA-Interim anomaly forcing is actually loaded (was previously a silent no-op)" begin
    tmpdir_anom = mktempdir()
    try
        clim_dir = joinpath(tmpdir_anom, "climatology")
        mkpath(clim_dir)
        write_field(name, value) = GREBClimate.jldopen(joinpath(clim_dir, name), "w") do file
            file["data"] = fill(value, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
            file["dim_names"] = ["lon", "lat", "time"]
        end

        # :rcp85 → CMIP5 RCP8.5 ensemble-mean anomaly
        write_field("cmip5.tsurf.rcp85.ensmean.forcing.jld2", 2.0)
        write_field("cmip5.zonal.wind.rcp85.ensmean.forcing.jld2", 3.0)
        write_field("cmip5.meridional.wind.rcp85.ensmean.forcing.jld2", 4.0)
        write_field("cmip5.omega.rcp85.ensmean.forcing.jld2", 5.0)
        write_field("cmip5.windspeed.rcp85.ensmean.forcing.jld2", 6.0)

        # :elnino / :lanina → ERA-Interim composite-mean anomaly
        for suffix in ("elnino", "lanina")
            write_field("erainterim.tsurf.$suffix.forcing.jld2", 7.0)
            write_field("erainterim.zonal.wind.$suffix.forcing.jld2", 8.0)
            write_field("erainterim.meridional.wind.$suffix.forcing.jld2", 9.0)
            write_field("erainterim.omega.$suffix.forcing.jld2", 10.0)
            write_field("erainterim.windspeed.$suffix.forcing.jld2", 11.0)
        end

        cfg = create_experiment_config(:rcp85)
        @test cfg.log_tsurf_ext && cfg.log_hwind_ext && cfg.log_omega_ext
        fields = ClimateFields()
        load_cc_anomaly_jld2!(tmpdir_anom, fields, cfg)
        @test all(==(2.0), fields.Tclim_anom_cc)
        @test all(==(3.0), fields.uclim_anom_cc)
        @test all(==(4.0), fields.vclim_anom_cc)
        @test all(==(5.0), fields.omegaclim_anom_cc)
        @test all(==(6.0), fields.wsclim_anom_cc)

        # init_model! applies the anomaly on top of the (here all-zero)
        # base climatology - Tclim must reflect it, not stay at zero.
        quiet() do
            init_model!(cfg, fields)
        end
        @test all(==(2.0), fields.Tclim)

        for (sym, suffix) in ((:elnino, "elnino"), (:lanina, "lanina"))
            cfg2 = create_experiment_config(sym)
            @test cfg2.log_tsurf_ext && cfg2.log_hwind_ext && cfg2.log_omega_ext
            fields2 = ClimateFields()
            load_enso_anomaly_jld2!(tmpdir_anom, fields2, cfg2, sym)
            @test all(==(7.0), fields2.Tclim_anom_enso)
            @test all(==(8.0), fields2.uclim_anom_enso)
            @test all(==(9.0), fields2.vclim_anom_enso)
            @test all(==(10.0), fields2.omegaclim_anom_enso)
            @test all(==(11.0), fields2.wsclim_anom_enso)
        end

        # Per-field gating: switching a gate off must not touch that
        # field even if its file is missing (no error, stays zero).
        cfg_partial = PhysicsConfig(experiment = :rcp85, log_tsurf_ext = true,
            log_hwind_ext = false, log_omega_ext = false)
        fields_partial = ClimateFields()
        load_cc_anomaly_jld2!(tmpdir_anom, fields_partial, cfg_partial)
        @test all(==(2.0), fields_partial.Tclim_anom_cc)
        @test all(==(0.0), fields_partial.uclim_anom_cc)
        @test all(==(0.0), fields_partial.omegaclim_anom_cc)

        # A missing required file must error loudly, not silently zero.
        rm(joinpath(clim_dir, "cmip5.tsurf.rcp85.ensmean.forcing.jld2"))
        @test_throws ErrorException load_cc_anomaly_jld2!(tmpdir_anom, ClimateFields(), cfg)
    finally
        rm(tmpdir_anom; recursive = true, force = true)
    end
end
