# Shortwave add-ons: the solar series and the stratospheric aerosol kernel.

@testset "load_solar_series returns multipliers relative to the table mean or a reference" begin
    with_tempdir() do dir
        path = joinpath(dir, "tsi.txt")
        write(path, "# year TSI
1990 1360.0

1991 1362.0
1992 1361.0
")
        s = load_solar_series(path)
        @test s[1991] ≈ Float32(1362 / 1361)
        @test s[1990] ≈ Float32(1360 / 1361)
        s2 = load_solar_series(path; reference = 1365.0)
        @test s2[1991] ≈ Float32(1362 / 1365)

        bad = joinpath(dir, "bad.txt")
        write(bad, "1990
")
        @test_throws ErrorException load_solar_series(bad)
        @test_throws ErrorException load_solar_series(joinpath(dir, "missing.txt"))
        empty_path = joinpath(dir, "empty.txt")
        write(empty_path, "# nothing
")
        @test_throws ErrorException load_solar_series(empty_path)
        @test_throws ErrorException load_solar_series(path; reference = -1.0)

        gap = joinpath(dir, "gap.txt")
        write(gap, "1990 1360.0
1992 1361.0
")
        @test_throws "gap" load_solar_series(gap)
        dup = joinpath(dir, "dup.txt")
        write(dup, "1990 1360.0
1990 1361.0
")
        @test_throws "more than once" load_solar_series(dup)
    end
end

@testset "apply_shortwave_addons! is an exact no-op unconfigured and applies the solar table otherwise" begin
    cfg = create_experiment_config(:full_model)
    mult = ones(Float32, Y)
    GREBClimate.apply_shortwave_addons!(mult, cfg, 1950, 1)
    @test all(==(1.0f0), mult)

    cfg.solar_scenario = Dict(1950 => 0.998f0)
    GREBClimate.apply_shortwave_addons!(mult, cfg, 1950, 1)
    @test all(==(0.998f0), mult)
    @test_throws ErrorException GREBClimate.apply_shortwave_addons!(mult, cfg, 1951, 1)
    @test_throws ErrorException GREBClimate.apply_shortwave_addons!(mult, cfg, 1950.5, 1)
end
