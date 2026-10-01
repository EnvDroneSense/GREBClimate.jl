# Grid constants, RunSpec, and resolve(): the hydrology coefficients and the
# tables a scenario reads.

@testset "grid constants" begin
    @test GREBClimate.xdim == 96
    @test GREBClimate.ydim == 48
    @test GREBClimate.nstep_yr == 730
end

@testset "RunSpec" begin
    @test RunSpec() == RunSpec(ctrl = 1, scnr = 1)
    @test fieldnames(RunSpec) == (:ctrl, :scnr)
end

@testset "resolve(::Hydrology): rain coefficients per scheme" begin
    coefficients(h) = (r = resolve(h); (r.c_q, r.c_rq, r.c_omega, r.c_omegastd))
    @test coefficients(Hydrology(rain = :original)) == (1.0f0, 0.0f0, 0.0f0, 0.0f0)
    @test coefficients(Hydrology(rain = :rh)) == (-1.391649f0, 3.018774f0, 0.0f0, 0.0f0)
    @test coefficients(Hydrology(rain = :omega)) == (0.862162f0, 0.0f0, -29.02096f0, 0.0f0)
    @test coefficients(Hydrology(rain = :rh_omega)) == (-0.2685845f0, 1.4591853f0, -26.9858807f0, 0.0f0)
    @test coefficients(Hydrology()) == (-1.88f0, 2.25f0, -17.69f0, 59.07f0)
    @test coefficients(Hydrology(rain_fit = :ncep)) == (-1.27f0, 1.99f0, -16.54f0, 21.15f0)
    h = resolve(Hydrology(rain = :rh, evaporation = :skin_gust))
    @test (h.rain, h.evaporation) == (:rh, :skin_gust)
end

@testset "resolve(::Config): tables come from the scenario" begin
    r = resolve(preset(:full_model))
    @test r.config == preset(:full_model)
    @test isempty(r.co2_table) && r.solar_table === nothing
    @test_throws ArgumentError resolve(preset(:custom_co2))          # no path
    with_tempdir() do dir
        write_ipcc_scenarios(dir, Dict("rcp45" => Dict(1950 => 401.0)))
        @test resolve(preset(:rcp45); jld2_dir = dir).co2_table == Dict(1950 => 401.0f0)
        path = joinpath(dir, "co2.txt")
        write(path, "1950 300\n")
        @test resolve(preset(:custom_co2; path)).co2_table == Dict(1950 => 300.0f0)
        write_solar_scenarios(dir)
        t = resolve(preset(:obliquity; index = 0); jld2_dir = dir).solar_table
        @test size(t) == (Y, N) && all(==(999.0f0), t)
    end
end
