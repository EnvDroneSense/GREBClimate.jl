# Processes, Hydrology, Corrections: validation, and the translation to the legacy
# PhysicsConfig. The legacy-comparison testsets go when PhysicsConfig does.

@testset "Processes, Hydrology, Corrections: defaults and validation" begin
    p = Processes()
    @test (p.atmosphere, p.clouds, p.hydrology, p.ocean, p.topography) == (true, :observed, :full, :full, :observed)
    @test Processes(clouds = :uniform).clouds === :uniform
    @test_throws ArgumentError Processes(clouds = :off)
    @test_throws ArgumentError Processes(ocean = :slab)
    @test_throws ArgumentError Hydrology(rain = :best)
    @test_throws ArgumentError Hydrology(rain = :original, rain_fit = :ncep)   # a fit only applies to :fitted
    @test mscm_hydrology() == Hydrology(rain = :original, evaporation = :original)
    @test Processes().co2 === true
    @test SpinUp(3).years == 3
    @test_throws ArgumentError SpinUp(-1)
end

@testset "Hydrology translates to the legacy log_rain/log_clim/log_eva codes" begin
    rains = (:original => -1, :fitted => 0, :rh => 1, :omega => 2, :rh_omega => 3)
    evas = (:original => -1, :skin => 0, :original_gust => 1, :skin_gust => 2)
    for (r, lr) in rains, (e, le) in evas, fit in (r === :fitted ? (:era, :ncep) : (:era,))
        cfg = GREBClimate._lower_hydrology!(PhysicsConfig(), Hydrology(rain = r, rain_fit = fit, evaporation = e))
        @test (cfg.log_rain, cfg.log_eva, cfg.log_clim) == (lr, le, fit === :ncep ? 1 : 0)
    end
end

# What the legacy switches actually do, read from the kernels and init_model!.
# Two configs with the same effect run the same model.
function switch_effect(c::PhysicsConfig)
    circulation = c.log_atmos_dmc && c.log_crcl_dmc && c.log_crcl_drsp
    return (atmosphere = c.log_atmos_dmc,
            circulation = circulation,
            transport_parts = circulation ? (c.log_hdif, c.log_hadv, c.log_vdif, c.log_vadv, c.log_conv) : nothing,
            hydro_kernel = c.log_atmos_dmc && c.log_hydro_dmc && c.log_hydro_drsp,
            humidity_update = c.log_hydro_dmc,
            cldclim = !c.log_clouds_drsp ? :uniform : !c.log_clouds_dmc ? :zero : :observed,
            qclim = !c.log_humid_drsp ? :uniform : !c.log_hydro_dmc ? :zero : :observed,
            ocean_capacity = c.log_ocean_dmc,          # cap_surf over ocean, seaice!
            mldclim = c.log_ocean_dmc ? (c.log_ocean_drsp ? :observed : :fixed) : nothing,
            deep_ocean = c.log_ocean_dmc && c.log_ocean_drsp,
            ice_albedo = c.log_ice,
            topography = c.log_topo_drsp ? :observed : :flat,
            corrections = GREBClimate._corrections_mode(c),
            co2 = c.log_co2_dmc)
end

# Legacy flags -> v2 options (test-only; the package only translates downwards)
function raise(c::PhysicsConfig)
    e = switch_effect(c)
    p = Processes(atmosphere = c.log_atmos_dmc,
                  clouds = e.cldclim === :zero ? :none : e.cldclim,
                  humidity = c.log_humid_drsp ? :observed : :uniform,
                  hydrology = !c.log_hydro_dmc ? :none : !c.log_hydro_drsp ? :no_evap_rain : :full,
                  ocean = !c.log_ocean_dmc ? :none : !c.log_ocean_drsp ? :mixed_layer : :full,
                  topography = e.topography, co2 = c.log_co2_dmc, ice_albedo = c.log_ice,
                  transport = c.log_crcl_dmc && c.log_crcl_drsp,
                  heat_diffusion = c.log_hdif, heat_advection = c.log_hadv,
                  vapour_diffusion = c.log_vdif, vapour_advection = c.log_vadv,
                  moisture_convergence = c.log_conv)
    corr = e.corrections === :stored ? Stored() : e.corrections === :spinup ? SpinUp(3) : NoCorrections()
    return p, corr
end

@testset "every legacy switch combination survives the translation ($name)" for (name, exp, keys) in (
        ("mean climate, 2048", :decon_mean_climate, (:log_clouds_dmc, :log_ocean_dmc, :log_atmos_dmc, :log_co2_dmc,
            :log_hydro_dmc, :log_qflux_dmc, :log_ice, :log_hdif, :log_hadv, :log_vdif, :log_vadv)),
        ("2xCO2 response, 1024", :decon_2xco2, (:log_topo_drsp, :log_clouds_drsp, :log_humid_drsp, :log_ocean_drsp,
            :log_hydro_drsp, :log_ice, :log_hdif, :log_hadv, :log_vdif, :log_vadv)))
    mismatches = 0
    for bits in Iterators.product(ntuple(_ -> (true, false), length(keys))...)
        legacy = create_experiment_config(exp; NamedTuple{keys}(bits)...)
        p, corr = raise(legacy)
        lowered = GREBClimate._lower_switches!(PhysicsConfig(experiment = exp), p, corr)
        mismatches += switch_effect(lowered) != switch_effect(legacy)
    end
    @test mismatches == 0
end

@testset "every topography and corrections combination translates" begin
    mode(p, c) = GREBClimate._corrections_mode(GREBClimate._lower_switches!(PhysicsConfig(), p, c))
    for topography in (:observed, :flat)
        p = Processes(topography = topography)
        @test mode(p, SpinUp(3)) === :spinup
        @test mode(p, Stored()) === :stored
        @test mode(p, NoCorrections()) === :none
    end
end
