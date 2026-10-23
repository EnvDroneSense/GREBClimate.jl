# greb_model! integration: runs end to end, what a run applies and when, and the
# flux corrections of the spin-up.

@testset "greb_model! baseline: default config runs to completion with the right output shape" begin
    result = quiet() do
        greb_model!(RunSpec(scnr = 0), preset(:full_model); jld2_dir = "", allow_uninitialized = true)
    end
    @test length(result.ctrl) == 12
    @test length(result.scnr) == 0
    @test result.ctrl[1] isa MonthlyRecord
end

@testset "qflux_correction! pulls Ts/To/q to climatology; Ta gets no correction" begin
    # All ocean, with a Ts climatology that varies from cell to cell
    fields = constant_fields(z_topo = -1.0)
    fields.cap_surf .= cap_ocean * 50.0f0
    fields.Ts_clim .= [280.0 + 5.0 * sin(i / 10.0) * cos(j / 8.0) for i in 1:X, j in 1:Y]
    fields.To_clim .= 279.0
    # Rain coefficients (1, 0, 0, 0): this hand-made climatology has no omega
    r = resolve(preset(:full_model; hydrology = (rain = :original,)))
    state = ModelState()
    ws = ModelWorkspace()

    Ts = fill(290.0, X, Y)
    Ta = fill(290.0, X, Y)
    q = fill(0.010, X, Y)
    To = fill(285.0, X, Y)

    qflux_correction!(340.0, Ts, Ta, q, To, fields, state, r, ws, 1)

    for name in (:Ts_flux_correction, :To_flux_correction, :q_flux_correction)
        @test any(!=(0.0), getfield(fields, name))
        @test all(isfinite, getfield(fields, name))
    end

    @test Ts ≈ fields.Ts_clim[:, :, 1]
    @test To ≈ fields.To_clim[:, :, 1]
    @test q ≈ fields.q_clim[:, :, 1]

    @test all(isfinite, Ta)
    @test Ta != fill(290.0, X, Y)
end

@testset "greb_model! swaps sw_solar for paleo experiments, restores after" begin
    # A paleo run's swapped solar table must not leak into a later run.
    fields = ClimateFields()
    saved_sw_solar = copy(fields.sw_solar)
    with_tempdir() do dir
        write_solar_scenarios(dir; sentinel = 999.0)

        cfg = preset(:paleo_231kyr; corrections = NoCorrections())
        @test all(==(999.0), resolve(cfg; jld2_dir = dir).solar_table)
        sw_solar(phase, run) = at_first_step(v -> copy(v.fields.sw_solar), run, cfg; phase, jld2_dir = dir, fields)
        # The control keeps the modern table, the scenario runs on the paleo one
        @test sw_solar(control, RunSpec(ctrl = 1, scnr = 0)) == saved_sw_solar
        @test all(==(999.0), sw_solar(scenario, RunSpec(ctrl = 0, scnr = 1)))

        # sw_solar restored to its pre-run value after greb_model! returns
        @test fields.sw_solar == saved_sw_solar
    end
end

@testset "a latitude CO2 mask applies in the scenario only" begin
    cfg = preset(:regional_co2_nh; corrections = NoCorrections())
    part(phase, run) = at_first_step(v -> copy(v.fields.co2_part), run, cfg; phase, jld2_dir = "")
    # The control runs on the full CO2 everywhere, as in the original code
    @test all(isone, part(control, RunSpec(ctrl = 1, scnr = 0)))
    scnr = part(scenario, RunSpec(ctrl = 0, scnr = 1))
    @test all(==(0.5f0), scnr[:, 1:24]) && all(isone, scnr[:, 25:48])
end

@testset "what a run logs: spin-up years are labelled, an empty scenario is not announced" begin
    cfg = preset(:full_model; corrections = SpinUp(1))
    logs(run) = Test.collect_test_logs() do
        greb_model!(run, cfg; fields = synthetic_fields(), allow_uninitialized = true)
    end[1]
    messages = [l.message for l in logs(RunSpec(ctrl = 1, scnr = 0))]
    @test count(startswith("spin-up year 1: Ts global mean"), messages) == 1
    @test count(startswith("1970: Ts global mean"), messages) == 1
    @test any(startswith("Control run"), messages) && !any(startswith("Scenario run"), messages)
end

@testset "a scenario runs on the CO2 or solar table its preset loads" begin
    with_tempdir() do dir
        write_ipcc_scenarios(dir, Dict("rcp45" => Dict(1950 => 401.0)))
        @test at_first_step(v -> v.CO2, RunSpec(ctrl = 0, scnr = 1),
                            preset(:rcp45; corrections = NoCorrections()); jld2_dir = dir) == 401.0f0

        co2_path = joinpath(dir, "my_co2.txt")
        write(co2_path, "1950 300.0\n")
        cfg = preset(:custom_co2; path = co2_path, corrections = NoCorrections())
        @test at_first_step(v -> v.CO2, RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = "") == 300.0f0
        # A missing path must raise a clear error, not silently default.
        @test_throws ArgumentError greb_model!(RunSpec(ctrl = 0, scnr = 1),
            preset(:custom_co2); jld2_dir = "", allow_uninitialized = true)

        write_solar_scenarios(dir)
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
    function run_mask(sym)
        f = ClimateFields()
        f.z_topo[1:(X - 48), :] .= 100.0f0   # left half land, right half ocean
        cfg = preset(sym; corrections = NoCorrections())
        result = run_synthetic(RunSpec(ctrl = 1, scnr = 0), cfg; fields = f)
        expected = ClimateFields()
        expected.z_topo .= f.z_topo
        apply_surface_mask!(cfg.scenario.co2_mask, expected, ice_climatology(result.ctrl))
        return (got = copy(f.co2_part), expected = expected.co2_part)
    end
    ocean = run_mask(:regional_co2_ocean)
    @test ocean.got == ocean.expected
    @test all(==(0.5f0), ocean.got[1:(X - 48), :])   # land halved whatever the ice
    land_ice = run_mask(:regional_co2_land_ice)
    @test land_ice.got == land_ice.expected
    @test all(isone, land_ice.got[1:(X - 48), :])    # land kept whatever the ice

    # Without a control run there is no ice cover to build the mask from
    @test_throws ArgumentError greb_model!(RunSpec(ctrl = 0, scnr = 1),
        preset(:regional_co2_ocean; corrections = NoCorrections()); jld2_dir = "", allow_uninitialized = true)
end

@testset "boundary anomalies: the add step adds the loaded anomaly, and no other" begin
    with_tempdir() do dir
        write_boundary_anomalies(dir)
        fields = ClimateFields()
        load_boundary_anomaly!(dir, fields, :cmip5_rcp85)
        # Adding an anomaly that is not the loaded one is an error, not a silent mix-up
        @test_throws ErrorException GREBClimate._add_boundary_anomaly!(BoundaryAnomaly(:elnino), fields)
        @test_throws ErrorException GREBClimate._add_boundary_anomaly!(BoundaryAnomaly(:elnino), ClimateFields())

        # The anomaly is added on top of the (here all-zero) base climatology
        quiet() do
            GREBClimate._add_boundary_anomaly!(BoundaryAnomaly(:cmip5_rcp85), fields)
        end
        @test all(==(2.0), fields.Ts_clim)

        for sym in (:elnino, :lanina)
            fields2 = ClimateFields()
            load_boundary_anomaly!(dir, fields2, sym)
            # Both composite files already carry their sign (the La Nina one is
            # a cold anomaly), so both experiments add them, as the Fortran does.
            quiet() do
                GREBClimate._add_boundary_anomaly!(BoundaryAnomaly(sym), fields2)
            end
            @test all(==(7.0), fields2.Ts_clim)
            @test all(==(8.0), fields2.u_clim)
        end
    end
end

@testset "boundary anomalies: a run adds them at the scenario start; the control runs on the plain climatology" begin
    with_tempdir() do dir
        write_boundary_anomalies(dir)
        for (p, value) in ((:rcp85_boundary, 2.0f0), (:lanina, 7.0f0))
            cfg = preset(p; corrections = NoCorrections())
            tclim(phase, run) = at_first_step(v -> copy(v.fields.Ts_clim), run, cfg; phase,
                                              jld2_dir = dir, fields = ClimateFields())
            @test all(iszero, tclim(control, RunSpec(ctrl = 1, scnr = 0)))
            @test all(==(value), tclim(scenario, RunSpec(ctrl = 0, scnr = 1)))
        end
    end
end

@testset "boundary anomalies: a second run reuses the loaded arrays; another source or directory reloads" begin
    # A second run on the same fields, directory and source does not read the
    # files again; another source or another directory does, into the same arrays.
    with_tempdir() do dir
        write_boundary_anomalies(dir)
        reused = ClimateFields()
        first_tclim(p, from) = at_first_step(v -> copy(v.fields.Ts_clim), RunSpec(ctrl = 0, scnr = 1),
                                            preset(p; corrections = NoCorrections()); phase = scenario,
                                            jld2_dir = from, fields = reused)
        @test all(==(7.0f0), first_tclim(:lanina, dir))
        @test reused.boundary_anomaly_source == (dir, :lanina)
        arrays = reused.boundary_anomaly
        reused.boundary_anomaly.Ts_clim .= 70.0f0              # a mark a reload would erase
        @test all(==(70.0f0), first_tclim(:lanina, dir))
        @test all(==(7.0f0), first_tclim(:elnino, dir))
        @test reused.boundary_anomaly_source == (dir, :elnino)
        @test all(==(2.0f0), first_tclim(:rcp85_boundary, dir))
        @test reused.boundary_anomaly === arrays               # another source reuses the five arrays
        reused.boundary_anomaly.Ts_clim .= 20.0f0
        @test all(==(20.0f0), first_tclim(:rcp85_boundary, dir))
        @test_throws ErrorException first_tclim(:rcp85_boundary, joinpath(dir, "elsewhere"))
        @test reused.boundary_anomaly_source == ("", :none)    # a failed load leaves no source
    end
end

@testset "flux-correction round trip: a run step from the spin-up's start lands on Ts_clim" begin
    # The spin-up solves for Ts_flux_correction so its step ends on Ts_clim; the run step
    # adds the same surface flux sum plus that Ts_flux_correction. A flux term added to
    # one of the two update loops but not the other moves Ts off Ts_clim here.
    f = synthetic_fields()
    r = resolve(preset(:full_model))
    ini = quiet(() -> init_model!(r, f))
    cap0 = copy(f.cap_surf)   # seaice! changes it during the spin-up year
    s = initial_state(ini)
    quiet() do
        qflux_correction!(ini.CO2_ctrl, s.Ts, s.Ta, s.q, s.To, f, ModelState(), r, ModelWorkspace(), 1)
    end
    function run_step()
        f.cap_surf .= cap0
        s = initial_state(ini)
        step!(control_time(), ini.CO2_ctrl, s, f, r)
        return s.Ts
    end
    @test isapprox(run_step(), f.Ts_clim[:, :, 1]; rtol = 1e-6)
    # Negative control: without the correction the same step does not land
    f.Ts_flux_correction[:, :, 1] .= 0.0f0
    @test !isapprox(run_step(), f.Ts_clim[:, :, 1]; rtol = 1e-6)
end

@testset "greb_model! restores the fields it changes, so one fields serves several runs" begin
    f = synthetic_fields()
    f.Ts_flux_correction .= 1.0f0   # stands in for loaded corrections; the run zeroes them here
    # sw_solar is not in this list: no part of this config swaps the solar
    # table, so the paleo testset is where its restore is checked.
    names = (:Ts_flux_correction, :q_flux_correction, :To_flux_correction, :z_topo, :cloud_clim, :q_clim, :mld_clim)
    before = Dict(n => copy(getfield(f, n)) for n in names)
    # Running without corrections zeroes them; flat topography, uniform clouds
    # and humidity and the mixed-layer ocean replace their climatologies.
    cfg = preset(:full_model; corrections = NoCorrections(),
                 processes = (topography = :flat, clouds = :uniform, humidity = :uniform, ocean = :mixed_layer))
    run_synthetic(RunSpec(ctrl = 1, scnr = 0), cfg; fields = f)
    for n in names
        @test getfield(f, n) == before[n]
    end
end

@testset "greb_model! spins up for as long as SpinUp says" begin
    # Ts after the first control step
    first_ts(years) = at_first_step(v -> copy(v.Ts), RunSpec(ctrl = 1, scnr = 0),
                                    preset(:full_model; corrections = SpinUp(years));
                                    phase = control, jld2_dir = "", fields = synthetic_fields())
    none, one_year = first_ts(0), first_ts(1)
    @test all(isfinite, none) && none != one_year
end

@testset "greb_model! takes the corrections it is given, with either topography" begin
    with_tempdir() do dir
        write_flux_corrections(dir; Ts = 0.5f0)
        for topography in (:observed, :flat)
            # Ts after the first control step; `preloaded` stands in for corrections already in `fields`
            function first_ts(c; preloaded = 0.0f0, jld2_dir = dir)
                f = synthetic_fields()
                f.Ts_flux_correction .= preloaded
                cfg = preset(:full_model; processes = (topography = topography,), corrections = c)
                return at_first_step(v -> copy(v.Ts), RunSpec(ctrl = 1, scnr = 0), cfg; phase = control, jld2_dir, fields = f)
            end
            stored = first_ts(Stored())
            @test all(isfinite, stored)
            # Stored reads the file: the same as a run on those values without a spin-up
            @test isequal(stored, first_ts(SpinUp(0); preloaded = 0.5f0, jld2_dir = ""))
            # Without a directory Stored keeps the corrections already in `fields`
            @test isequal(stored, first_ts(Stored(); preloaded = 0.5f0, jld2_dir = ""))
            # NoCorrections zeroes whatever was there
            @test isequal(first_ts(NoCorrections(); preloaded = 0.5f0), first_ts(SpinUp(0); jld2_dir = ""))
            @test !isequal(stored, first_ts(NoCorrections()))
            # SpinUp computes them
            @test !isequal(stored, first_ts(SpinUp(1)))
        end
    end
    # A directory without the corrections file is an error, not a run on zeros
    with_tempdir() do empty_dir
        @test_throws ArgumentError run_synthetic(RunSpec(ctrl = 1, scnr = 0),
                                                 preset(:full_model; corrections = Stored()); jld2_dir = empty_dir)
    end
end

@testset "greb_model! reloads Stored corrections only when fields does not hold them from the file" begin
    with_tempdir() do dir
        write_flux_corrections(dir; Ts = 0.5f0)
        cfg = preset(:full_model; corrections = Stored())
        first_ts(f) = at_first_step(v -> copy(v.Ts), RunSpec(ctrl = 1, scnr = 0), cfg;
                                    phase = control, jld2_dir = dir, fields = f)

        file = joinpath(dir, "climatology", "flux_corrections.jld2")

        # A run on fields without the marker loads the file for the run and leaves
        # the arrays and the marker as it found them
        f = synthetic_fields()
        fresh = first_ts(f)
        @test f.flux_source == ""
        @test all(iszero, f.Ts_flux_correction)

        # Loaded explicitly, the marker is set and a run keeps the arrays
        quiet() do
            load_flux_corrections!(dir, f)
        end
        @test f.flux_source == file
        @test isequal(first_ts(f), fresh)
        @test f.flux_source == file
        @test all(==(0.5f0), f.Ts_flux_correction)

        # Arrays changed under a set marker are used as they are: no reload
        f.Ts_flux_correction .= 0.7f0
        @test !isequal(first_ts(f), fresh)
        @test all(==(0.7f0), f.Ts_flux_correction)

        # Resetting the marker reloads the file for the run, then restores the arrays
        f.flux_source = ""
        @test isequal(first_ts(f), fresh)
        @test all(==(0.7f0), f.Ts_flux_correction)
    end
end

@testset "greb_model! returns the time of every record" begin
    result = run_synthetic(RunSpec(ctrl = 2, scnr = 2), Config(scenario = Scenario(start_year = 1850), corrections = NoCorrections()))
    @test length(result.ctrl) == length(result.ctrl_time) == 24
    @test length(result.scnr) == length(result.scnr_time) == 24
    @test result.ctrl_time == record_times(1970, 24)
    @test result.scnr_time == record_times(1850, 24)   # the scenario's start year

    @testset "scnr_anomaly says what scnr holds" begin
        @test result.scnr_anomaly                       # default output is :anomaly, both runs have records
        absolute = run_synthetic(RunSpec(ctrl = 1, scnr = 1),
                             Config(scenario = Scenario(output = :absolute), corrections = NoCorrections()))
        @test !absolute.scnr_anomaly
        @test !run_synthetic(RunSpec(ctrl = 0, scnr = 1), preset(:co2_double; corrections = NoCorrections())).scnr_anomaly
    end

    @testset "scenario without control" begin
        result = run_synthetic(RunSpec(ctrl = 0, scnr = 1), preset(:co2_double; corrections = NoCorrections()))
        @test isempty(result.ctrl) && isempty(result.ctrl_time)
        @test result.scnr_time == record_times(1950, 12)
    end

    @testset "years at and below zero" begin
        config = Config(scenario = Scenario(start_year = -1, output = :absolute), corrections = NoCorrections())
        result = run_synthetic(RunSpec(ctrl = 0, scnr = 2), config)
        @test result.scnr_time == record_times(-1, 24)
    end
end
