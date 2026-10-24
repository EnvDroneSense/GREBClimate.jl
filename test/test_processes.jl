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
        ini = quiet(() -> init_model!(f, resolve(config)))
        return f, ini
    end
    f, ini = init((;))
    @test f.z_topo[1, 1] == 500.0f0 && all(==(1.0f0), f.Ts_flux_correction) && ini.CO2_ctrl == 340
    @test all(==(0.0f0), init((clouds = :none,))[1].cloud_clim)
    @test all(==(0.7f0), init((clouds = :uniform,))[1].cloud_clim)
    @test all(==(0.0f0), init((hydrology = :none,))[1].q_clim)
    @test all(==(0.0052f0), init((hydrology = :none, humidity = :uniform))[1].q_clim)   # uniform wins
    @test all(==(d_ocean), init((ocean = :mixed_layer,))[1].mld_clim)
    @test all(==(cap_land), init((ocean = :none,))[1].cap_surf)
    @test init((topography = :flat,))[1].z_topo[1, 1] == 1.0f0
    @test init((co2 = false,))[2].CO2_ctrl == 0
    @test all(==(0.0f0), init((;); corrections = NoCorrections())[1].Ts_flux_correction)
end

@testset "init_model! changes only fields that greb_model! restores or derives" begin
    # A write that the override table and `_mutated_fields` do not know would leak into a reused `fields`
    derived = (:dTrad, :z_ocean, :wz_air, :wz_vapor, :rain_limit, :cap_surf,
               :u_clim_pos, :u_clim_neg, :v_clim_pos, :v_clim_neg, :co2_part)
    options = ((;), (topography = :flat,), (clouds = :none,), (clouds = :uniform,), (hydrology = :none,),
               (humidity = :uniform,), (ocean = :mixed_layer,),
               (topography = :flat, clouds = :uniform, hydrology = :none, humidity = :uniform, ocean = :mixed_layer))
    for processes in options, corrections in (SpinUp(1), NoCorrections())
        config = Config(; processes = Processes(; processes...), corrections)
        f = synthetic_fields()
        f.z_topo[1, 1] = 500.0f0
        f.Ts_flux_correction .= 1.0f0
        before = Dict(n => copy(getfield(f, n)) for n in fieldnames(ClimateFields) if getfield(f, n) isa Array)
        quiet(() -> init_model!(f, resolve(config)))
        changed = [n for (n, a) in before if getfield(f, n) != a]
        @test issubset(changed, union(derived, _mutated_fields(config)))
    end
end
