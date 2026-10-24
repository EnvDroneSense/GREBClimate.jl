# JLD2 loading, dataset resolution, and converter/archive consistency.

@testset "load_climatology: every file is read; a missing flux-corrections file or table is an error" begin
    with_tempdir() do tmpdir
        map2(name, v) = write_field(joinpath(tmpdir, "static", name), v, (X, Y))
        clim(name, v) = write_field(joinpath(tmpdir, "climatology", name), v)
        solar(v, dims = (Y, N)) = write_field(joinpath(tmpdir, "solar", "solar_radiation.clim.jld2"), v, dims;
                                              dim_names = ["lat", "time"])
        map2("global.topography.jld2", 1.0)
        map2("greb.glaciers.jld2", 2.0)
        clim("ncep.tsurf.1948-2007.clim.jld2", 3.0)
        clim("ncep.zonal_wind.850hpa.clim.jld2", 4.0)
        clim("ncep.meridional_wind.850hpa.clim.jld2", 5.0)
        clim("ncep.atmospheric_humidity.clim.jld2", 6.0)
        clim("ncep.soil_moisture.clim.jld2", 7.0)
        clim("isccp.cloud_cover.clim.jld2", 8.0)
        clim("woce.ocean_mixed_layer_depth.clim.jld2", 9.0)
        clim("Tocean.clim.jld2", 10.0)
        clim("erainterim.omega.vertmean.clim.jld2", 11.0)
        clim("erainterim.omega_std.vertmean.clim.jld2", 12.0)
        clim("erainterim.windspeed.850hpa.clim.jld2", 13.0)
        solar(14.0)

        # No flux-corrections file yet: an error that names it, not zeros
        @test_throws "flux_corrections.jld2" load_climatology(tmpdir; dataset = :ncep)
        @test_throws ArgumentError load_flux_corrections!(ClimateFields(), tmpdir)
        # `corrections = false` loads the rest and leaves the corrections at zero
        without = load_climatology(tmpdir; dataset = :ncep, corrections = false)
        @test without.loaded && all(==(3.0), without.Ts_clim)
        @test all(iszero, without.Ts_flux_correction) && all(iszero, without.q_flux_correction) && all(iszero, without.To_flux_correction)

        # A dataset name the loader does not know is an error, not NCEP
        @test_throws ArgumentError load_climatology(tmpdir; dataset = :era5)

        # A file that lacks one of the three tables is an error too
        write_flux_corrections(tmpdir; Ts = 15.0, q = nothing, To = nothing)
        @test_throws "vapour_flux_correction" load_climatology(tmpdir; dataset = :ncep)

        write_flux_corrections(tmpdir; Ts = 15.0, q = 16.0, To = 17.0)

        fields = load_climatology(tmpdir; dataset = :ncep)
        @test all(==(1.0), fields.z_topo)
        @test all(==(2.0), fields.glacier)
        @test all(==(3.0), fields.Ts_clim)
        @test all(==(4.0), fields.u_clim)
        @test all(==(14.0), fields.sw_solar)
        @test all(==(15.0), fields.Ts_flux_correction)
        @test all(==(16.0), fields.q_flux_correction)
        @test all(==(17.0), fields.To_flux_correction)

        # Tasks reading the same file at once all get its content
        path = joinpath(tmpdir, "climatology", "Tocean.clim.jld2")
        reads = [Threads.@spawn read_field(path).data for _ in 1:8]
        @test all(t -> all(==(10.0), fetch(t)), reads)

        # A solar table of the wrong shape is reported, with the file's name
        solar(0.0, (Y, 2))
        @test_throws "solar_radiation.clim.jld2" load_climatology(tmpdir; dataset = :ncep)
    end

    with_tempdir() do missing_parent
        @test_throws ArgumentError load_climatology(joinpath(missing_parent, "nonexistent"))
    end
end

@testset "IPCC scenario CO2 tables load under the right on-disk key" begin
    # The resolve step loads each table once. Assert the loader against every
    # preset's key instead of paying a simulated year each.
    expected = Dict(:rcp26 => 400.0, :rcp45 => 401.0, :rcp60 => 402.0,
                    :ssp119 => 300.0, :ssp126 => 301.0, :ssp245 => 302.0,
                    :ssp460 => 303.0, :ssp585 => 304.0, :historical_co2 => 280.73)
    key(p) = preset(p).scenario.co2.key
    on_disk(p) = p === :rcp60 ? "rcp6" : string(key(p))   # the dataset's name for RCP6.0
    with_tempdir() do dir
        write_ipcc_scenarios(dir, Dict(on_disk(p) => Dict(1950 => co2) for (p, co2) in expected))
        for (p, co2) in expected
            @test isapprox(load_co2_scenario(dir, key(p))[1950], co2; atol = 1e-3)
        end

        # The table's name is :rcp60; the dataset's own key is not accepted
        @test key(:rcp60) === :rcp60
        @test_throws ArgumentError load_co2_scenario(dir, :rcp6)
    end
end

@testset "custom CO2 trajectory loads from a plain-text file" begin
    with_tempdir() do dir
        co2_path = joinpath(dir, "my_co2.txt")
        write(co2_path, "# comment line, should be skipped\n1950 300.0\n1951 301.0\n\n")

        # The parser is the thing under test - assert it directly.
        @test load_co2_custom(co2_path) == Dict(1950 => 300.0, 1951 => 301.0)

        # A malformed line must raise without leaving the file open: on
        # Windows an open handle makes the rm below fail.
        bad_path = joinpath(dir, "bad_co2.txt")
        write(bad_path, "1950 300.0\n1951\n")
        @test_throws ArgumentError load_co2_custom(bad_path)
        rm(bad_path)
        @test !isfile(bad_path)
    end
end

@testset "paleo/orbital solar tables load" begin
    with_tempdir() do dir
        write_solar_scenarios(dir)
        # The loader is the mechanism; assert all three tables directly.
        for kind in (:paleo, :obliquity, :eccentricity)
            table = load_solar_forcing(dir, kind, 0)
            @test size(table) == (Y, N)
            @test all(==(999.0f0), table)
        end
        # An orbital index that is not in the table
        @test_throws ArgumentError load_solar_forcing(dir, :obliquity, 7)
        @test_throws ArgumentError load_solar_forcing(dir, :eccentricity, 7)
        @test_throws ArgumentError load_solar_forcing(dir, :bogus)
    end
end

@testset "boundary anomalies: each source's files load into BoundaryAnomalyFields; a missing file is an error" begin
    with_tempdir() do dir
        expected = write_boundary_anomalies(dir)
        fields = ClimateFields()
        @test fields.boundary_anomaly === nothing      # no anomaly arrays until one is loaded
        load_boundary_anomaly!(fields, dir, :cmip5_rcp85)
        @test fields.boundary_anomaly isa BoundaryAnomalyFields
        @test fields.boundary_anomaly_source == (dir, :cmip5_rcp85)
        for (source, want) in expected
            loaded = ClimateFields()
            load_boundary_anomaly!(loaded, dir, source)
            for name in fieldnames(BoundaryAnomalyFields)
                @test all(==(want[name]), getfield(loaded.boundary_anomaly, name))
            end
        end

        # A missing required file must error loudly, not silently zero.
        rm(joinpath(dir, "climatology", "cmip5.tsurf.rcp85.ensmean.forcing.jld2"))
        @test_throws ArgumentError load_boundary_anomaly!(ClimateFields(), dir, :cmip5_rcp85)
        @test_throws ArgumentError load_boundary_anomaly!(ClimateFields(), dir, :neutral)
    end
end

@testset "greb_data_dir resolution order" begin
    with_tempdir() do tmp_a
        with_tempdir() do tmp_b
            withenv("GREB_DATA" => tmp_b) do
                # explicit path wins over everything
                @test greb_data_dir(tmp_a) == tmp_a
                # ...and over the environment even with allow_download off
                @test greb_data_dir(tmp_a; allow_download = false) == tmp_a
                # GREB_DATA wins over the repo-local dataset
                @test greb_data_dir() == tmp_b
            end

            withenv("GREB_DATA" => nothing) do
                # a non-existent explicit path is an error, not a silent fallback
                @test_throws ErrorException greb_data_dir(joinpath(tmp_a, "nope"))
            end
            # so is a GREB_DATA pointing nowhere
            withenv("GREB_DATA" => joinpath(tmp_a, "nope")) do
                @test_throws ErrorException greb_data_dir()
            end

            # A downloaded dataset is the directory `<load path>/GREB-input-data`.
            # Stand one up so the cache step runs on machines that never downloaded it.
            # DataDeps splits DATADEPS_LOAD_PATH on ':', which cuts a Windows drive
            # letter off, so the load path is given relative to the working directory.
            cd(tmp_b) do
                withenv("GREB_DATA" => nothing, "DATADEPS_LOAD_PATH" => ".", "DATADEPS_NO_STANDARD_LOAD_PATH" => "true") do
                    @test GREBClimate._cached_datadep_path() === nothing
                    cache = mkpath(joinpath(pwd(), GREBClimate.DATA_DEP_NAME))
                    @test abspath(GREBClimate._cached_datadep_path()) == cache
                    # the repo-local dataset wins over the cache; without it the cache is used
                    local_dir = normpath(joinpath(@__DIR__, "..", "greb_input_data"))
                    @test abspath(greb_data_dir(; allow_download = false)) == (isdir(local_dir) ? local_dir : cache)
                end
            end

            withenv("GREB_DATA" => nothing) do
                @test greb_data_dir(""; allow_download = false) == greb_data_dir(; allow_download = false)
            end
        end
    end
end

@testset "published dataset archive constants are coherent" begin
    @test occursin(r"^[0-9a-f]{64}$", GREBClimate.DATA_SHA256)

    data_src = read(joinpath(@__DIR__, "..", "src", "data.jl"), String)
    @test match(r"const DATA_RELEASE_TAG = \"([^\"]+)\"", data_src).captures[1] ==
          GREBClimate.DATA_RELEASE_TAG
end

@testset "dataset file list: 33 single-field files, all of them in the dataset tools' list" begin
    names = GREBClimate.dataset_field_files()
    @test length(names) == 33
    # Both ENSO events and both datasets are in it
    @test "erainterim.omega.lanina.forcing" in names && "erainterim.tsurf.elnino.forcing" in names
    @test "ncep.tsurf.1948-2007.clim" in names && "erainterim.tsurf.1979-2015.clim" in names
    # Every name fills an array that ClimateFields has
    lists = (GREBClimate._STATIC_FILES, values(GREBClimate._CLIMATOLOGY_FILES)..., GREBClimate._COMMON_CLIMATOLOGY_FILES,
             GREBClimate._FLUX_CORRECTION_KEYS)
    @test all(l -> all(f -> hasfield(ClimateFields, f), keys(l)), lists)
    # Every source of a boundary anomaly names one file per anomaly array
    @test all(s -> keys(GREBClimate._boundary_anomaly_files(s)) == fieldnames(BoundaryAnomalyFields),
              GREBClimate._BOUNDARY_ANOMALY_SOURCES)

    # The tools keep their own list: it may hold more, never less
    tools = Module()
    Base.include(tools, joinpath(@__DIR__, "..", "tools", "dataset", "fields.jl"))
    @test issubset(names, tools.MODEL_FIELD_NAMES)
    @test issubset(GREBClimate._COMBINED_FILES, tools.COMBINED_FILE_NAMES)
    @test Set(tools.FLUX_CORRECTION_NAMES) == Set(values(GREBClimate._FLUX_CORRECTION_KEYS))
end
