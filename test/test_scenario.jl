# Config, Scenario and preset(): every preset translates to exactly the legacy
# PhysicsConfig of the experiment it replaces. The legacy comparison goes when
# PhysicsConfig does; test_presets.jl keeps guarding the forcing after that.

# new preset name => legacy experiment
const RENAMED = Dict(:rcp85_boundary => :rcp85, :co2_abrupt_reverse => :co2_step, :a1b => :a1b_scenario)
legacy_name(p) = get(RENAMED, p, p)

# Corrections compare by what the run does: a lowered config names them, a
# legacy one leaves them to its switches
function same_legacy(a::PhysicsConfig, b::PhysicsConfig; skip = ())
    names = [n for n in fieldnames(PhysicsConfig) if !(n in skip) && n !== :flux_corrections]
    bad = [n for n in names if getfield(a, n) != getfield(b, n)]
    GREBClimate._corrections_mode(a) === GREBClimate._corrections_mode(b) || push!(bad, :flux_corrections)
    isempty(bad) || @info "differing PhysicsConfig fields" bad
    return isempty(bad)
end

@testset "preset names: renames applied, :constant_topo dropped" begin
    names = preset_names()
    @test :rcp85_boundary in names && :co2_abrupt_reverse in names && :a1b in names
    @test !(:rcp85 in names) && !(:co2_step in names) && !(:a1b_scenario in names) && !(:constant_topo in names)
    @test length(names) == 41     # 42 legacy presets, minus :constant_topo
    @test_throws ArgumentError preset(:constant_topo)
    @test_throws ArgumentError preset(:rcp85)
    @test_throws ArgumentError preset(:not_a_preset)
end

@testset "preset(:$p) lowers to the legacy config" for p in filter(p -> !(p in (:custom_co2, :decon_mean_climate, :decon_2xco2)), preset_names())
    @test same_legacy(GREBClimate._lower(preset(p)), create_experiment_config(legacy_name(p)))
end

@testset "parameterised presets carry their parameter" begin
    @test same_legacy(GREBClimate._lower(preset(:eccentricity; index = 32)),
                      create_experiment_config(:eccentricity; orbital_index = 32))
    @test same_legacy(GREBClimate._lower(preset(:obliquity; index = 95)),
                      create_experiment_config(:obliquity; orbital_index = 95))
    @test same_legacy(GREBClimate._lower(preset(:earth_sun_distance; pct = 2.5)),
                      create_experiment_config(:earth_sun_distance; earth_sun_distance_pct = 2.5))
    @test same_legacy(GREBClimate._lower(preset(:custom_co2; path = "co2.txt")),
                      create_experiment_config(:custom_co2; co2_path = "co2.txt"))
end

@testset "deconstruction presets" begin
    # The legacy :decon_* symbols run the same forcing as :full_model / :co2_double.
    # The mean-climate preset runs the MSCM physics on the stored corrections.
    mscm_legacy(; kw...) = (c = create_experiment_config(:decon_mean_climate; kw...);
                            c.log_rain = -1; c.log_eva = -1; c.log_conv = false; c.flux_corrections = :stored; c)
    skip = (:experiment,)
    @test same_legacy(GREBClimate._lower(preset(:decon_mean_climate)), mscm_legacy(); skip)
    @test same_legacy(GREBClimate._lower(preset(:decon_mean_climate; processes = (ocean = :none, clouds = :none))),
                      mscm_legacy(log_ocean_dmc = false, log_clouds_dmc = false); skip)
    @test same_legacy(GREBClimate._lower(preset(:decon_2xco2; processes = Processes(topography = :flat, humidity = :uniform),
                                                corrections = Stored())),
                      create_experiment_config(:decon_2xco2; log_topo_drsp = false, log_humid_drsp = false); skip)
end

@testset "preset defaults: corrections, and changes to the preset's physics" begin
    c = preset(:decon_mean_climate)
    @test (c.processes, c.hydrology, c.corrections) == (Processes(moisture_convergence = false), mscm_hydrology(), Stored())
    # A NamedTuple changes the preset's physics; a Processes/Hydrology replaces it
    @test preset(:decon_mean_climate; processes = (ocean = :none,)).processes ==
          Processes(moisture_convergence = false, ocean = :none)
    @test preset(:decon_mean_climate; hydrology = (rain = :fitted,)).hydrology == Hydrology(evaporation = :original)
    @test preset(:decon_mean_climate; processes = Processes(ocean = :none)).processes == Processes(ocean = :none)
    @test preset(:co2_double; processes = (topography = :flat,)).processes == Processes(topography = :flat)
    # Flat topography spins up like everything else unless told otherwise
    @test preset(:co2_double; processes = (topography = :flat,)).corrections == SpinUp(3)
    @test preset(:decon_2xco2).corrections == SpinUp(3)
    @test_throws ArgumentError preset(:co2_double; processes = (oceans = :none,))
end

@testset "CO2 switched off: 0 ppm in the control and the scenario" begin
    fields = ClimateFields()
    ice = zeros(Float32, X, Y, 12)
    for p in (:decon_mean_climate, :full_model, :co2_double, :rcp45, :sst_plus1)
        cfg = GREBClimate._lower(preset(p; processes = (co2 = false,)))
        @test quiet(() -> init_model!(cfg, fields)).CO2_ctrl == 0
        cfg.co2_scenario = Dict(1950 => 400f0)
        @test forcing(1, 1950, cfg, fields, ice).CO2 == 0
    end
    # The solar part of the forcing is untouched
    cfg = GREBClimate._lower(preset(:solar_plus27; processes = (co2 = false,)))
    @test forcing(1, 1950, cfg, fields, ice).sw_solar_forcing ≈ 1392 / 1365
end

@testset "scenarios: one per legacy experiment; others are refused for now" begin
    s = [preset(p).scenario for p in preset_names() if !(p in (:decon_mean_climate, :decon_2xco2))]
    @test allunique(s)
    @test_throws ArgumentError GREBClimate._lower(Config(scenario = Scenario(co2 = ConstantCO2(500))))
    @test_throws ArgumentError GREBClimate._lower(Config(scenario = Scenario(co2 = ConstantCO2(680), solar = SolarConstant(27))))
    @test Config().scenario == preset(:full_model).scenario
    @test Config() == preset(:full_model)
end
