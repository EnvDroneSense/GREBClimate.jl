# Processes, Hydrology, Corrections: validation, and what each Processes
# option replaces in init_model!.

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

@testset "init_model!: what each Processes option replaces" begin
    function init(processes; corrections = SpinUp(3))
        f = synthetic_fields()
        f.z_topo[1, 1] = 500.0f0
        f.Ts_flux_correction .= 1.0f0
        config = Config(; processes = Processes(; processes...), corrections)
        ini = quiet(() -> init_model!(resolve(config), f))
        return f, ini
    end
    f, ini = init((;))
    @test f.z_topo[1, 1] == 500.0f0 && all(==(1.0f0), f.Ts_flux_correction) && ini.CO2_ctrl == 340
    @test all(==(0.0f0), init((clouds = :none,))[1].cloud_clim)
    @test all(==(0.7f0), init((clouds = :uniform,))[1].cloud_clim)
    @test all(==(0.0f0), init((hydrology = :none,))[1].q_clim)
    @test all(==(0.0052f0), init((hydrology = :none, humidity = :uniform))[1].q_clim)   # uniform wins
    @test all(==(GREBClimate.d_ocean), init((ocean = :mixed_layer,))[1].mld_clim)
    @test all(==(GREBClimate.cap_land), init((ocean = :none,))[1].cap_surf)
    @test init((topography = :flat,))[1].z_topo[1, 1] == 1.0f0
    @test init((co2 = false,))[2].CO2_ctrl == 0
    @test all(==(0.0f0), init((;); corrections = NoCorrections())[1].Ts_flux_correction)
end
