# Config, Scenario and preset(): names, defaults, parameters, and forcing from
# free combinations of scenario parts. test_presets.jl guards what every preset
# imposes against the stored reference.

# A config resolved with a stand-in CO2 table, without the dataset
resolved(c; table=Dict(1950 => 400.0f0)) = ResolvedConfig(c, resolve(c.hydrology), table, nothing)

@testset "preset() refuses a name that is not a preset" begin
    # The set of names itself is pinned by the preset reference (test_presets.jl)
    @test_throws ArgumentError preset(:not_a_preset)
end

@testset "parameterised presets carry their parameter" begin
    # The default rows (32, 95) are pinned by the preset reference. Use other
    # rows here, so an ignored `index` shows
    @test SolarTable(:paleo).index == 0
    @test preset(:eccentricity; index=5).scenario.solar == SolarTable(:eccentricity, 5)
    @test preset(:obliquity; index=10).scenario.solar == SolarTable(:obliquity, 10)
    @test preset(:earth_sun_distance; percent=2.5).scenario.solar == EarthSunDistance(2.5)
    @test preset(:custom_co2; path="co2.txt").scenario.co2 == CO2File("co2.txt")
end

@testset "preset defaults: corrections, and changes to the preset's physics" begin
    c = preset(:decon_mean_climate)
    @test (c.processes, c.hydrology, c.corrections) == (Processes(moisture_convergence=false), mscm_hydrology(), Stored())
    # A NamedTuple changes the preset's physics; a Processes/Hydrology replaces it
    @test preset(:decon_mean_climate; processes=(ocean=:none,)).processes ==
          Processes(moisture_convergence=false, ocean=:none)
    @test preset(:decon_mean_climate; hydrology=(rain=:fitted,)).hydrology == Hydrology(evaporation=:original)
    @test preset(:decon_mean_climate; processes=Processes(ocean=:none)).processes == Processes(ocean=:none)
    @test preset(:co2_double; processes=(topography=:flat,)).processes == Processes(topography=:flat)
    # Flat topography spins up like everything else unless told otherwise
    @test preset(:co2_double; processes=(topography=:flat,)).corrections == SpinUp(3)
    # The response deconstruction runs the MSCM physics too, on its own spin-up
    d = preset(:decon_2xco2)
    @test (d.processes, d.hydrology, d.corrections) == (Processes(moisture_convergence=false), mscm_hydrology(), SpinUp(3))
    @test_throws ArgumentError preset(:co2_double; processes=(oceans=:none,))
end

@testset "CO2 switched off: 0 ppm in the control and the scenario" begin
    fields = ClimateFields()
    for p in (:decon_mean_climate, :full_model, :co2_double, :rcp45, :sst_plus1)
        r = resolved(preset(p; processes=(co2=false,)))
        @test quiet(() -> init_model!(fields, r)).CO2_ctrl == 0
        @test forcing(scenario_time(1950, 1), r).CO2 == 0
    end
    # The solar part of the forcing is untouched
    @test forcing(scenario_time(1950, 1), resolved(preset(:solar_plus27; processes=(co2=false,)))).sw_solar_forcing ≈ 1392 / 1365
end

@testset "free combinations of scenario parts" begin
    @test allunique(preset(p).scenario for p in preset_names() if !(p in (:decon_mean_climate, :decon_2xco2)))
    @test Config() == preset(:full_model)
    # Parts no preset combines
    r = resolved(Config(scenario=Scenario(co2=ConstantCO2(500), solar=SolarConstant(27))))
    @test forcing(scenario_time(1950, 1), r) == (CO2=500.0f0, sw_solar_forcing=(1365.0f0 + 27.0f0) / 1365.0f0)
    r = resolved(Config(scenario=Scenario(co2=CO2Step(300, 900, 2000), solar=SolarCycle(2, 22))))
    @test forcing(scenario_time(1999, 1), r).CO2 == 300 && forcing(scenario_time(2000, 1), r).CO2 == 900
    @test forcing(scenario_time(2000, 1), r).sw_solar_forcing ≈ (1365 + 2sin(2π * 2000 / 22)) / 1365
    r = resolved(Config(scenario=Scenario(co2=SeasonalCO2(1000, 300, :boreal_summer))))
    @test forcing(scenario_time(1950, 1), r).CO2 == 300 && forcing(scenario_time(1950, 300), r).CO2 == 1000
    @test_throws ErrorException forcing(scenario_time(1951, 1), resolved(preset(:rcp45)))    # year not in the table
end

@testset "A1B ramp: 700 ppm in 2100, held afterwards" begin
    r = resolved(preset(:a1b))
    @test forcing(scenario_time(1950, 1), r).CO2 == 310 && forcing(scenario_time(2050, 1), r).CO2 == 520
    @test forcing(scenario_time(2100, 1), r).CO2 == 700
    @test forcing(scenario_time(2101, 1), r).CO2 == 700 && forcing(scenario_time(2300, 1), r).CO2 == 700
end

@testset "SeasonalCO2 season is half the year, from 1 October" begin
    r = resolve(Config(scenario=Scenario(co2=SeasonalCO2(680, 340, :boreal_winter))))
    co2(step) = forcing(scenario_time(1950, step), r).CO2
    @test GREBClimate._winter_first_step == first_step_of(10, 1) == 547
    @test GREBClimate._winter_last_step == first_step_of(4, 1) == 181
    @test co2.([1, 181, 182, 546, 547, 730]) == Float32[680, 680, 340, 340, 680, 680]
    @test count(==(680.0f0), co2.(1:nstep_yr)) == nstep_yr ÷ 2
    @test co2(731) == 680.0f0     # the second year starts in winter again
end

@testset "CO2 masks: latitude bands, surface masks from the annual-mean ice cover" begin
    f = ClimateFields()
    apply_co2_mask!(f, LatitudeMask(:nh))
    @test all(==(0.5f0), f.co2_part[:, 1:24]) && all(isone, f.co2_part[:, 25:48])
    apply_co2_mask!(f, UniformMask())          # resets a reused fields
    @test all(isone, f.co2_part)
    # Each band keeps the full CO2 at the latitudes it names (longitude 1 is
    # clear of the edge rows' every-fourth-longitude exception)
    lat = lat_grid
    full(band) = (apply_co2_mask!(f, LatitudeMask(band)); lat[f.co2_part[1, :] .== 1.0f0])
    @test full(:nh) == lat[lat .> 0]
    @test full(:sh) == lat[lat .< 0]
    @test full(:tropics) == lat[-33.75 .< lat .< 30]
    @test full(:extratropics) == lat[(lat .< -33.75) .| (lat .> 30)]
    # Every fourth longitude of the two halved rows at the band edge keeps
    # the full CO2
    apply_co2_mask!(f, LatitudeMask(:tropics))
    @test findall(isone, f.co2_part[:, 15]) == 4:4:X && findall(isone, f.co2_part[:, 33]) == 4:4:X
    apply_co2_mask!(f, LatitudeMask(:extratropics))
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

    apply_surface_mask!(fields, SurfaceMask(:ocean), icmn_ctrl)
    @test fields.co2_part[1, 1] == 1.0  # January said "ice"; annual mean says no
    @test fields.co2_part[2, 1] == 0.5  # January said "no ice"; annual mean says yes

    # The land/ice variant inverts the ocean mask and exempts ice cells.
    f_li = ClimateFields()
    apply_surface_mask!(f_li, SurfaceMask(:land_ice), icmn_ctrl)
    @test f_li.co2_part[1, 1] == 0.5  # ocean cell, not annual-mean ice
    @test f_li.co2_part[2, 1] == 1.0  # annual-mean ice -> exempted back to 1.0

    # Every other mask is a no-op here.
    for mask in (UniformMask(), LatitudeMask(:sh))
        f_noop = ClimateFields()
        apply_surface_mask!(f_noop, mask, icmn_ctrl)
        @test all(isone, f_noop.co2_part)
    end
end
