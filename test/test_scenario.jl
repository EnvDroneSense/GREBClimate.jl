# Config, Scenario and preset(): names, defaults, parameters, and forcing from
# free combinations of scenario parts. test_presets.jl guards what every preset
# imposes against the stored reference.

# A config resolved with a stand-in CO2 table, without the dataset
resolved(c; table = Dict(1950 => 400.0f0)) = ResolvedConfig(c, resolve(c.hydrology), table, nothing)

@testset "preset() refuses a name that is not a preset" begin
    # The set of names itself is pinned by the preset reference (test_presets.jl)
    @test_throws ArgumentError preset(:not_a_preset)
end

@testset "orbital presets default to the near-modern table row" begin
    @test preset(:eccentricity).scenario.solar == SolarTable(:eccentricity, 32)
    @test preset(:obliquity).scenario.solar == SolarTable(:obliquity, 95)
    @test SolarTable(:obliquity) == SolarTable(:obliquity, 95)
    @test SolarTable(:paleo).index == 0
end

@testset "parameterised presets carry their parameter" begin
    # not the default rows (32, 95): an ignored `index` must show
    @test preset(:eccentricity; index = 5).scenario.solar == SolarTable(:eccentricity, 5)
    @test preset(:obliquity; index = 10).scenario.solar == SolarTable(:obliquity, 10)
    @test preset(:earth_sun_distance; pct = 2.5).scenario.solar == EarthSunDistance(2.5)
    @test preset(:custom_co2; path = "co2.txt").scenario.co2 == CO2File("co2.txt")
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
    # The response deconstruction runs the MSCM physics too, on its own spin-up
    d = preset(:decon_2xco2)
    @test (d.processes, d.hydrology, d.corrections) == (Processes(moisture_convergence = false), mscm_hydrology(), SpinUp(3))
    @test_throws ArgumentError preset(:co2_double; processes = (oceans = :none,))
end

@testset "CO2 switched off: 0 ppm in the control and the scenario" begin
    fields = ClimateFields()
    for p in (:decon_mean_climate, :full_model, :co2_double, :rcp45, :sst_plus1)
        r = resolved(preset(p; processes = (co2 = false,)))
        @test quiet(() -> init_model!(r, fields)).CO2_ctrl == 0
        @test forcing(1, 1950, r).CO2 == 0
    end
    # The solar part of the forcing is untouched
    @test forcing(1, 1950, resolved(preset(:solar_plus27; processes = (co2 = false,)))).sw_solar_forcing ≈ 1392 / 1365
end

@testset "free combinations of scenario parts" begin
    @test allunique(preset(p).scenario for p in preset_names() if !(p in (:decon_mean_climate, :decon_2xco2)))
    @test Config() == preset(:full_model)
    # Parts no preset combines
    r = resolved(Config(scenario = Scenario(co2 = ConstantCO2(500), solar = SolarConstant(27))))
    @test forcing(1, 1950, r) == (CO2 = 500.0f0, sw_solar_forcing = (1365.0f0 + 27.0f0) / 1365.0f0)
    r = resolved(Config(scenario = Scenario(co2 = CO2Step(300, 900, 2000), solar = SolarCycle(2, 22))))
    @test forcing(1, 1999, r).CO2 == 300 && forcing(1, 2000, r).CO2 == 900
    @test forcing(1, 2000, r).sw_solar_forcing ≈ (1365 + 2sin(2π * 2000 / 22)) / 1365
    r = resolved(Config(scenario = Scenario(co2 = SeasonalCO2(1000, 300, :boreal_summer))))
    @test forcing(1, 1950, r).CO2 == 300 && forcing(300, 1950, r).CO2 == 1000
    @test_throws ErrorException forcing(1, 1951, resolved(preset(:rcp45)))    # year not in the table
end

@testset "SeasonalCO2 season is half the year, from 1 October" begin
    r = resolve(Config(scenario=Scenario(co2=SeasonalCO2(680, 340, :boreal_winter))))
    co2(step) = GREBClimate.forcing(step, 1950, r).CO2
    @test GREBClimate._winter_first_step == GREBClimate.first_step_of(10, 1) == 547
    @test GREBClimate._winter_last_step == GREBClimate.first_step_of(4, 1) == 181
    @test co2.([1, 181, 182, 546, 547, 730]) == Float32[680, 680, 340, 340, 680, 680]
    @test count(==(680.0f0), co2.(1:nstep_yr)) == nstep_yr ÷ 2
    @test co2(731) == 680.0f0     # the second year starts in winter again
end
