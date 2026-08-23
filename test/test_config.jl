# Grid constants, PhysicsConfig, create_experiment_config, hydrology parameters.

@testset "grid constants" begin
    @test GREBClimate.xdim == 96
    @test GREBClimate.ydim == 48
    @test GREBClimate.nstep_yr == 730
end

@testset "PhysicsConfig" begin
    for exp in (:full_model, :constant_topo, :co2_double, :co2_quadruple,
                :elnino, :lanina, :rcp26, :rcp45, :rcp60, :rcp85, :ssp119,
                :ssp126, :ssp245, :ssp460, :ssp585, :historical_co2,
                :decon_mean_climate, :decon_2xco2)
        c = create_experiment_config(exp)
        @test c.experiment == exp
    end

    @test create_experiment_config(:co2_double).co2_concentration == 680.0
end

@testset "create_experiment_config: :custom_co2 and decon presets (§7.2/§7.3)" begin
    cfg_custom = create_experiment_config(:custom_co2; co2_path="/tmp/my_co2.txt")
    @test cfg_custom.experiment == :custom_co2
    @test cfg_custom.custom_co2_path == "/tmp/my_co2.txt"
    @test create_experiment_config(:custom_co2).custom_co2_path == ""

    # decon_mean_climate: defaults all true, one override propagates
    cfg_dmc = create_experiment_config(:decon_mean_climate)
    for switch in (cfg_dmc.log_clouds_dmc, cfg_dmc.log_ocean_dmc, cfg_dmc.log_atmos_dmc,
                   cfg_dmc.log_co2_dmc, cfg_dmc.log_hydro_dmc, cfg_dmc.log_qflux_dmc,
                   cfg_dmc.log_ice, cfg_dmc.log_hdif, cfg_dmc.log_hadv,
                   cfg_dmc.log_vdif, cfg_dmc.log_vadv)
        @test switch == true
    end
    cfg_dmc_off = create_experiment_config(:decon_mean_climate; log_ocean_dmc=false)
    @test cfg_dmc_off.log_ocean_dmc == false
    @test cfg_dmc_off.log_clouds_dmc == true  # untouched switches stay at default

    # decon_2xco2: defaults all true + doubled CO2, one override propagates
    cfg_drsp = create_experiment_config(:decon_2xco2)
    @test cfg_drsp.co2_concentration == 680.0
    for switch in (cfg_drsp.log_topo_drsp, cfg_drsp.log_clouds_drsp, cfg_drsp.log_humid_drsp,
                   cfg_drsp.log_ocean_drsp, cfg_drsp.log_hydro_drsp,
                   cfg_drsp.log_ice, cfg_drsp.log_hdif, cfg_drsp.log_hadv,
                   cfg_drsp.log_vdif, cfg_drsp.log_vadv)
        @test switch == true
    end
    cfg_drsp_off = create_experiment_config(:decon_2xco2; log_topo_drsp=false)
    @test cfg_drsp_off.log_topo_drsp == false
    @test cfg_drsp_off.log_clouds_drsp == true  # untouched switches stay at default
end

@testset "create_experiment_config covers every experiment the model dispatches on" begin
    all_experiments = (
        :full_model, :constant_topo, :a1b_scenario,
        :co2_double, :co2_quadruple, :co2_10x, :co2_half, :co2_zero,
        :co2_sine_wave, :co2_step, :solar_plus27, :solar_cycle_11yr,
        :paleo_231kyr, :paleo_solar_modern_co2, :modern_solar_paleo_co2,
        :obliquity, :eccentricity, :earth_sun_distance,
        :elnino, :lanina, :rcp26, :rcp45, :rcp60, :rcp85,
        :ssp119, :ssp126, :ssp245, :ssp460, :ssp585,
        :historical_co2, :custom_co2, :sst_plus1,
        :regional_co2_nh, :regional_co2_sh, :regional_co2_tropics,
        :regional_co2_extratropics, :regional_co2_ocean,
        :regional_co2_land_ice, :regional_co2_winter, :regional_co2_summer,
        :decon_mean_climate, :decon_2xco2,
    )
    @test length(all_experiments) == 42
    for exp in all_experiments
        cfg = create_experiment_config(exp)
        @test cfg.experiment === exp
    end
    @test_throws ErrorException create_experiment_config(:no_such_experiment)

    for exp in (:co2_10x, :co2_half, :co2_zero, :co2_sine_wave, :co2_step,
                :a1b_scenario, :solar_cycle_11yr, :sst_plus1,
                :regional_co2_nh, :regional_co2_ocean, :regional_co2_winter,
                :paleo_solar_modern_co2, :modern_solar_paleo_co2)
        factory = create_experiment_config(exp)
        manual = PhysicsConfig(experiment = exp)
        for f in fieldnames(PhysicsConfig)
            @test getfield(factory, f) == getfield(manual, f)
        end
    end

    # Orbital / distance experiments plumb their parameter through.
    @test create_experiment_config(:obliquity; orbital_index = 3).orbital_index == 3
    @test create_experiment_config(:eccentricity; orbital_index = 7).orbital_index == 7
    @test create_experiment_config(:earth_sun_distance;
              earth_sun_distance_pct = 1.5).earth_sun_distance_pct == 1.5f0

    # Unchanged presets stay unchanged.
    @test create_experiment_config(:co2_double).co2_concentration == 680.0f0
    @test create_experiment_config(:co2_quadruple).co2_concentration == 1360.0f0
    @test create_experiment_config(:paleo_231kyr).co2_concentration == 200.0f0
    @test create_experiment_config(:constant_topo).log_topo_drsp == false
    @test create_experiment_config(:elnino).log_tsurf_ext == true
    @test create_experiment_config(:rcp85).log_omega_ext == true

    # A log_* keyword aimed at a branch that ignores it now warns instead
    # of vanishing silently.
    @test_logs (:warn,) create_experiment_config(:co2_double; log_ice = false)
    @test create_experiment_config(:co2_double; log_ice = false).log_ice == true
end

@testset "set_hydrology_parameters! matches HYDRO_PARAMS for every log_rain value" begin
    for log_rain in (-1, 0, 1, 2, 3)
        cfg = create_experiment_config(:full_model)
        cfg.log_rain = log_rain
        set_hydrology_parameters!(cfg)
        expected = GREBClimate.HYDRO_PARAMS[log_rain]
        @test (cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd) == expected
    end

    cfg = create_experiment_config(:full_model)
    cfg.log_rain = 0
    cfg.log_clim = 1
    set_hydrology_parameters!(cfg)
    @test all(isapprox.((cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd), (-1.27, 1.99, -16.54, 21.15); atol = 1e-4))
end

@testset "set_hydrology_parameters! errors on invalid log_rain" begin
    # An out-of-range log_rain must fail loudly, not silently fall back
    # to a default parameter set.
    cfg = create_experiment_config(:full_model)
    cfg.log_rain = 4
    @test_throws ErrorException set_hydrology_parameters!(cfg)
end
