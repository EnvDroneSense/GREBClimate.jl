# greb_model! integration: presets end to end, scenario tables, CO2 masks, flux corrections.

@testset "greb_model! baseline: default config runs to completion with the right output shape" begin
    result = quiet() do
        greb_model!(RunSpec(scnr = 0), preset(:full_model); jld2_dir = "", allow_uninitialized = true)
    end
    @test length(result.ctrl) == 12
    @test length(result.scnr) == 0
    @test result.ctrl[1] isa MonthlyRecord
end

@testset "hydro! is finite under every evaporation and rain scheme" begin
    # Each evaporation scheme is a separate `@turbo` block in `hydro!`; the
    # default pair runs end to end in the baseline testset.
    fields = synthetic_fields()
    Ts = fill(288.0f0, X, Y)
    q = fill(0.005f0, X, Y)
    for (evaporation, rain) in ((:original, :original), (:skin, :rh), (:original_gust, :omega), (:skin_gust, :rh_omega))
        r = resolve(preset(:full_model; hydrology = (; evaporation, rain)))
        quiet() do
            init_model!(r, fields)
        end
        out = hydro!(Ts, q, fields, TimeState(1, 1), Processes(), r.hydrology, CirculationWorkspace())
        @test all(isfinite, out.Q_lat)
        @test all(isfinite, out.dq_rain)
    end
end

@testset "qflux_correction! pulls Ts/To/q to climatology; Ta gets no correction" begin
    # All ocean (topography 0 m), so it needs a mixed layer
    fields = ClimateFields()
    fields.mldclim .= 50.0
    fields.cap_surf .= GREBClimate.cap_ocean * 50.0f0
    for j in 1:GREBClimate.ydim, i in 1:GREBClimate.xdim
        fields.Tclim[i, j, :] .= 280.0 + 5.0 * sin(i / 10.0) * cos(j / 8.0)
        fields.Toclim[i, j, :] .= 279.0
        fields.qclim[i, j, :] .= 0.006
    end
    # Rain coefficients (1, 0, 0, 0): this hand-made climatology has no omega
    r = resolve(preset(:full_model; hydrology = (rain = :original,)))
    state = ModelState()
    ts = TimeState(1, 1)
    ws = CirculationWorkspace()

    Ts = fill(290.0, GREBClimate.xdim, GREBClimate.ydim)
    Ta = fill(290.0, GREBClimate.xdim, GREBClimate.ydim)
    q = fill(0.010, GREBClimate.xdim, GREBClimate.ydim)
    To = fill(285.0, GREBClimate.xdim, GREBClimate.ydim)

    GREBClimate.qflux_correction!(340.0, Ts, Ta, q, To, fields, state, ts, r, ws, 1)

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

        cfg = preset(:paleo_231kyr; corrections = NoCorrections())
        @test all(==(distinctive_value), resolve(cfg; jld2_dir = tmpdir).solar_table)
        sw_solar(phase, run) = at_first_step(v -> copy(v.fields.sw_solar), run, cfg; phase, jld2_dir = tmpdir, fields)
        # The control keeps the modern table, the scenario runs on the paleo one
        @test sw_solar(:ctrl, RunSpec(ctrl = 1, scnr = 0)) == saved_sw_solar
        @test all(==(distinctive_value), sw_solar(:scnr, RunSpec(ctrl = 0, scnr = 1)))

        # sw_solar restored to its pre-run value after greb_model! returns
        @test fields.sw_solar == saved_sw_solar
    finally
        rm(tmpdir; recursive = true, force = true)
    end
end

@testset "a latitude CO2 mask applies in the scenario only" begin
    cfg = preset(:regional_co2_nh; corrections = NoCorrections())
    part(phase, run) = at_first_step(v -> copy(v.fields.co2_part), run, cfg; phase, jld2_dir = "")
    # The control runs on the full CO2 everywhere, as in the original code
    @test all(isone, part(:ctrl, RunSpec(ctrl = 1, scnr = 0)))
    scnr = part(:scnr, RunSpec(ctrl = 0, scnr = 1))
    @test all(==(0.5f0), scnr[:, 1:24]) && all(isone, scnr[:, 25:48])
end

@testset "CO2 masks: latitude bands, surface masks from the annual-mean ice cover" begin
    f = ClimateFields()
    GREBClimate.apply_co2_mask!(LatitudeMask(:nh), f)
    @test all(==(0.5f0), f.co2_part[:, 1:24]) && all(isone, f.co2_part[:, 25:48])
    GREBClimate.apply_co2_mask!(UniformMask(), f)          # resets a reused fields
    @test all(isone, f.co2_part)
    # Each band keeps the full CO2 at the latitudes it names (longitude 1 is
    # clear of the edge rows' every-fourth-longitude exception)
    lat = GREBClimate.lat_grid
    full(band) = (GREBClimate.apply_co2_mask!(LatitudeMask(band), f); lat[f.co2_part[1, :] .== 1.0f0])
    @test full(:nh) == lat[lat .> 0]
    @test full(:sh) == lat[lat .< 0]
    @test full(:tropics) == lat[-33.75 .< lat .< 30]
    @test full(:extratropics) == lat[(lat .< -33.75) .| (lat .> 30)]
    # Every fourth longitude of the two halved rows at the band edge keeps
    # the full CO2
    GREBClimate.apply_co2_mask!(LatitudeMask(:tropics), f)
    @test findall(isone, f.co2_part[:, 15]) == 4:4:X && findall(isone, f.co2_part[:, 33]) == 4:4:X
    GREBClimate.apply_co2_mask!(LatitudeMask(:extratropics), f)
    @test findall(isone, f.co2_part[:, 16]) == 4:4:X && findall(isone, f.co2_part[:, 32]) == 4:4:X

    fields = ClimateFields()  # z_topo defaults to 0 everywhere -> land branch never fires
    icmn_ctrl = zeros(Float64, X, Y, 12)
    # Cell A: January alone >= 0.5, but the other 11 months are 0 ->
    # annual mean ~0.083, NOT ice under the Fortran-matching rule.
    icmn_ctrl[1, 1, 1] = 1.0
    # Cell B: January alone < 0.5, but the other 11 months are 1.0 ->
    # annual mean ~0.917, IS ice under the Fortran-matching rule.
    icmn_ctrl[2, 1, 1] = 0.0
    icmn_ctrl[2, 1, 2:12] .= 1.0

    GREBClimate.apply_dynamic_co2_mask!(SurfaceMask(:ocean), fields, icmn_ctrl)
    @test fields.co2_part[1, 1] == 1.0  # January said "ice"; annual mean says no
    @test fields.co2_part[2, 1] == 0.5  # January said "no ice"; annual mean says yes

    # The land/ice variant inverts the ocean mask and exempts ice cells.
    f_li = ClimateFields()
    GREBClimate.apply_dynamic_co2_mask!(SurfaceMask(:land_ice), f_li, icmn_ctrl)
    @test f_li.co2_part[1, 1] == 0.5  # ocean cell, not annual-mean ice
    @test f_li.co2_part[2, 1] == 1.0  # annual-mean ice -> exempted back to 1.0

    # Every other mask is a no-op here.
    for mask in (UniformMask(), LatitudeMask(:sh))
        f_noop = ClimateFields()
        GREBClimate.apply_dynamic_co2_mask!(mask, f_noop, icmn_ctrl)
        @test all(isone, f_noop.co2_part)
    end
end

@testset "IPCC scenario CO2 tables load under the right on-disk key" begin
    # The resolve step loads each table once. Assert the loader against every
    # preset's key instead of paying a simulated year each.
    expected = Dict(:rcp26 => 400.0, :rcp45 => 401.0, :rcp60 => 402.0,
                    :ssp119 => 300.0, :ssp126 => 301.0, :ssp245 => 302.0,
                    :ssp460 => 303.0, :ssp585 => 304.0, :historical_co2 => 280.73)
    key(p) = preset(p).scenario.co2.key
    with_tempdir() do dir
        write_ipcc_scenarios(dir, Dict(string(key(p)) => Dict(1950 => co2) for (p, co2) in expected))
        for (p, co2) in expected
            @test isapprox(load_co2_scenario_jld2(dir, key(p))[1950], co2; atol = 1e-3)
        end

        # The scenario runs on the table's value
        @test at_first_step(v -> v.CO2, RunSpec(ctrl = 0, scnr = 1),
                            preset(:rcp45; corrections = NoCorrections()); jld2_dir = dir) == Float32(expected[:rcp45])
    end
end

@testset "custom CO2 trajectory loads from a plain-text file" begin
    with_tempdir() do dir
        co2_path = joinpath(dir, "my_co2.txt")
        write(co2_path, "# comment line, should be skipped\n1950 300.0\n1951 301.0\n\n")

        # The parser is the thing under test - assert it directly.
        @test load_custom_co2_scenario(co2_path) == Dict(1950 => 300.0, 1951 => 301.0)

        # A malformed line must raise without leaving the file open: on
        # Windows an open handle makes the rm below fail.
        bad_path = joinpath(dir, "bad_co2.txt")
        write(bad_path, "1950 300.0\n1951\n")
        @test_throws ErrorException load_custom_co2_scenario(bad_path)
        rm(bad_path)
        @test !isfile(bad_path)

        # ...and the scenario runs on the file's value.
        cfg = preset(:custom_co2; path = co2_path, corrections = NoCorrections())
        @test at_first_step(v -> v.CO2, RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = "") == 300.0f0

        # A missing path must raise a clear error, not silently default.
        @test_throws ArgumentError greb_model!(RunSpec(ctrl = 0, scnr = 1),
            preset(:custom_co2); jld2_dir = "", allow_uninitialized = true)
    end
end

@testset "paleo/orbital solar tables load; an orbital scenario runs on its table" begin
    with_tempdir() do dir
        write_solar_scenarios(dir)
        # The loader is the mechanism; assert all three tables directly.
        for kind in (:paleo, :obliquity, :eccentricity)
            table = load_solar_forcing_jld2(dir, kind, 0)
            @test size(table) == (Y, N)
            @test all(==(999.0f0), table)
        end

        cfg = preset(:obliquity; index = 0, corrections = NoCorrections())
        @test all(==(999.0f0), at_first_step(v -> copy(v.fields.sw_solar), RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = dir))
    end
end

@testset "greb_model!: sst_plus1, the deconstruction presets, the dynamic regional mask" begin
    # sst_plus1: the ocean surface is held 1 K above the climatology (all
    # ocean and 0 K here) and the CO2 stays at the control value
    seen = at_first_step(v -> (Ts = copy(v.Ts), CO2 = v.CO2), RunSpec(ctrl = 0, scnr = 1),
                         preset(:sst_plus1; corrections = NoCorrections()); point = :after_tendencies, jld2_dir = "")
    @test all(==(1.0f0), seen.Ts) && seen.CO2 == 340

    # The original-GREB physics of the deconstruction presets, end to end
    result_dmc = quiet() do
        greb_model!(RunSpec(ctrl = 1, scnr = 0),
                    preset(:decon_mean_climate);
                    jld2_dir = "", allow_uninitialized = true)
    end
    @test length(result_dmc.ctrl) == 12
    @test at_first_step(v -> v.CO2, RunSpec(ctrl = 0, scnr = 1),
                        preset(:decon_2xco2; corrections = NoCorrections()); jld2_dir = "") == 680

    # greb_model! builds the mask from the control run's own ice cover. The
    # data-free control is not physical, so compare against the mask that ice
    # cover gives rather than against a fixed pattern.
    run_mask(sym) = begin
        f = ClimateFields()
        f.z_topo[1:(X - 48), :] .= 100.0f0   # left half land, right half ocean
        cfg = preset(sym; corrections = NoCorrections())
        result = quiet() do
            greb_model!(RunSpec(ctrl = 1, scnr = 0), cfg;
                        jld2_dir = "", fields = f, allow_uninitialized = true)
        end
        expected = ClimateFields()
        expected.z_topo .= f.z_topo
        GREBClimate.apply_dynamic_co2_mask!(cfg.scenario.co2_mask, expected, compute_annual_ice_climatology(result.ctrl))
        (got = copy(f.co2_part), expected = expected.co2_part)
    end
    ocean = run_mask(:regional_co2_ocean)
    @test ocean.got == ocean.expected
    @test all(==(0.5f0), ocean.got[1:(X - 48), :])   # land halved whatever the ice
    land_ice = run_mask(:regional_co2_land_ice)
    @test land_ice.got == land_ice.expected
    @test all(isone, land_ice.got[1:(X - 48), :])    # land kept whatever the ice
end

@testset "boundary anomalies: the files load and are added to the climatology, in the scenario only" begin
    tmpdir_anom = mktempdir()
    try
        clim_dir = joinpath(tmpdir_anom, "climatology")
        mkpath(clim_dir)
        write_field(name, value) = GREBClimate.jldopen(joinpath(clim_dir, name), "w") do file
            file["data"] = fill(value, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
            file["dim_names"] = ["lon", "lat", "time"]
        end

        # :rcp85_boundary: CMIP5 RCP8.5 ensemble-mean anomaly
        write_field("cmip5.tsurf.rcp85.ensmean.forcing.jld2", 2.0)
        write_field("cmip5.zonal.wind.rcp85.ensmean.forcing.jld2", 3.0)
        write_field("cmip5.meridional.wind.rcp85.ensmean.forcing.jld2", 4.0)
        write_field("cmip5.omega.rcp85.ensmean.forcing.jld2", 5.0)
        write_field("cmip5.windspeed.rcp85.ensmean.forcing.jld2", 6.0)

        # :elnino / :lanina: ERA-Interim composite-mean anomaly
        for suffix in ("elnino", "lanina")
            write_field("erainterim.tsurf.$suffix.forcing.jld2", 7.0)
            write_field("erainterim.zonal.wind.$suffix.forcing.jld2", 8.0)
            write_field("erainterim.meridional.wind.$suffix.forcing.jld2", 9.0)
            write_field("erainterim.omega.$suffix.forcing.jld2", 10.0)
            write_field("erainterim.windspeed.$suffix.forcing.jld2", 11.0)
        end

        fields = ClimateFields()
        load_cc_anomaly_jld2!(tmpdir_anom, fields)
        @test all(==(2.0), fields.Tclim_anom_cc)
        @test all(==(3.0), fields.uclim_anom_cc)
        @test all(==(4.0), fields.vclim_anom_cc)
        @test all(==(5.0), fields.omegaclim_anom_cc)
        @test all(==(6.0), fields.wsclim_anom_cc)

        # The scenario-start step applies the anomaly on top of the (here
        # all-zero) base climatology - Tclim must reflect it, not stay at zero.
        quiet() do
            GREBClimate._apply_boundary_anomalies!(BoundaryAnomaly(:cmip5_rcp85), fields)
        end
        @test all(==(2.0), fields.Tclim)

        for (sym, suffix) in ((:elnino, "elnino"), (:lanina, "lanina"))
            fields2 = ClimateFields()
            load_enso_anomaly_jld2!(tmpdir_anom, fields2, sym)
            @test all(==(7.0), fields2.Tclim_anom_enso)
            @test all(==(8.0), fields2.uclim_anom_enso)
            @test all(==(9.0), fields2.vclim_anom_enso)
            @test all(==(10.0), fields2.omegaclim_anom_enso)
            @test all(==(11.0), fields2.wsclim_anom_enso)

            # Both composite files already carry their sign (the La Nina one is
            # a cold anomaly), so both experiments add them, as the Fortran does.
            quiet() do
                GREBClimate._apply_boundary_anomalies!(BoundaryAnomaly(sym), fields2)
            end
            @test all(==(7.0), fields2.Tclim)
            @test all(==(8.0), fields2.uclim)
        end

        # In a run the anomaly reaches the climatology at the start of the
        # scenario; the control runs on the plain one.
        for (p, anomaly) in ((:rcp85_boundary, 2.0f0), (:lanina, 7.0f0))
            cfg = preset(p; corrections = NoCorrections())
            tclim(phase, run) = at_first_step(v -> copy(v.fields.Tclim), run, cfg; phase,
                                              jld2_dir = tmpdir_anom, fields = ClimateFields())
            @test all(iszero, tclim(:ctrl, RunSpec(ctrl = 1, scnr = 0)))
            @test all(==(anomaly), tclim(:scnr, RunSpec(ctrl = 0, scnr = 1)))
        end

        # A missing required file must error loudly, not silently zero.
        rm(joinpath(clim_dir, "cmip5.tsurf.rcp85.ensmean.forcing.jld2"))
        @test_throws ErrorException load_cc_anomaly_jld2!(tmpdir_anom, ClimateFields())
    finally
        rm(tmpdir_anom; recursive = true, force = true)
    end
end

@testset "flux-correction round trip: a run step from the spin-up's start lands on Tclim" begin
    # The spin-up solves for TF_correct so its step ends on Tclim; the run step
    # adds the same surface flux sum plus that TF_correct. A flux term added to
    # one of the two update loops but not the other moves Ts off Tclim here.
    f = synthetic_fields()
    r = resolve(preset(:full_model))
    ini = quiet(() -> init_model!(r, f))
    cap0 = copy(f.cap_surf)   # seaice! changes it during the spin-up year
    start() = (copy(ini.Ts_ini), copy(ini.Ta_ini), copy(ini.q_ini), copy(ini.To_ini))
    Ts, Ta, q, To = start()
    quiet() do
        qflux_correction!(ini.CO2_ctrl, Ts, Ta, q, To, f, ModelState(), TimeState(1, 1), r,
                          CirculationWorkspace(), 1)
    end
    function run_step()
        f.cap_surf .= cap0
        Ts, Ta, q, To = start()
        quiet() do
            time_loop!(1, 1970, ini.CO2_ctrl, 1, 0, Ts, Ta, q, To, MonthlyRecord[], f, ModelState(),
                       CirculationWorkspace(), MonthlyAccumulator(), TimeState(1, 1), r)
        end
        return Ts
    end
    @test isapprox(run_step(), f.Tclim[:, :, 1]; rtol = 1e-6)
    # Negative control: without the correction the same step does not land
    f.TF_correct[:, :, 1] .= 0.0f0
    @test !isapprox(run_step(), f.Tclim[:, :, 1]; rtol = 1e-6)
end

@testset "greb_model! restores the fields it changes, so one fields serves several runs" begin
    f = synthetic_fields()
    f.TF_correct .= 1.0f0   # stands in for loaded corrections; the run zeroes them here
    # sw_solar is not in this list: no part of this config swaps the solar
    # table, so the paleo testset is where its restore is checked.
    names = (:TF_correct, :qF_correct, :ToF_correct, :z_topo, :cldclim, :qclim, :mldclim)
    before = Dict(n => copy(getfield(f, n)) for n in names)
    # Running without corrections zeroes them; flat topography, uniform clouds
    # and humidity and the mixed-layer ocean replace their climatologies.
    cfg = preset(:full_model; corrections = NoCorrections(),
                 processes = (topography = :flat, clouds = :uniform, humidity = :uniform, ocean = :mixed_layer))
    quiet() do
        greb_model!(RunSpec(ctrl = 1, scnr = 0), cfg; jld2_dir = "", fields = f,
                    allow_uninitialized = true)
    end
    for n in names
        @test getfield(f, n) == before[n]
    end
end

@testset "greb_model! spins up for as long as SpinUp says" begin
    # Ts after the first control step
    first_ts(years) = at_first_step(v -> copy(v.Ts), RunSpec(ctrl = 1, scnr = 0),
                                    preset(:full_model; corrections = SpinUp(years));
                                    phase = :ctrl, jld2_dir = "", fields = synthetic_fields())
    none, one_year = first_ts(0), first_ts(1)
    @test all(isfinite, none) && none != one_year
end

@testset "greb_model! takes the corrections it is given, with either topography" begin
    with_tempdir() do dir
        mkpath(joinpath(dir, "climatology"))
        GREBClimate.jldopen(joinpath(dir, "climatology", "flux_corrections.jld2"), "w") do f
            f["Tsurf_flux_correction"] = fill(0.5f0, X, Y, GREBClimate.nstep_yr)
            f["vapour_flux_correction"] = zeros(Float32, X, Y, GREBClimate.nstep_yr)
            f["Tocean_flux_correction"] = zeros(Float32, X, Y, GREBClimate.nstep_yr)
        end
        for topography in (:observed, :flat)
            # Ts after the first control step; `preloaded` stands in for corrections already in `fields`
            function first_ts(c; preloaded = 0.0f0, jld2_dir = dir)
                f = synthetic_fields()
                f.TF_correct .= preloaded
                cfg = preset(:full_model; processes = (topography = topography,), corrections = c)
                return at_first_step(v -> copy(v.Ts), RunSpec(ctrl = 1, scnr = 0), cfg; phase = :ctrl, jld2_dir, fields = f)
            end
            stored = first_ts(Stored())
            @test all(isfinite, stored)
            # Stored reads the file: the same as a run on those values without a spin-up
            @test isequal(stored, first_ts(SpinUp(0); preloaded = 0.5f0, jld2_dir = ""))
            # NoCorrections zeroes whatever was there
            @test isequal(first_ts(NoCorrections(); preloaded = 0.5f0), first_ts(SpinUp(0); jld2_dir = ""))
            @test !isequal(stored, first_ts(NoCorrections()))
            # SpinUp computes them
            @test !isequal(stored, first_ts(SpinUp(1)))
        end
    end
end
