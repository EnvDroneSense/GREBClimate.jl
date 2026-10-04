# JLD2 loading, dataset resolution, and converter/archive consistency.

@testset "load_greb_jld2!: every file is read; a missing flux-corrections file or table is an error" begin
    write2(path, v) = (mkpath(dirname(path)); GREBClimate.jldopen(path, "w") do f
        f["data"] = fill(v, GREBClimate.xdim, GREBClimate.ydim); f["dim_names"] = ["lon", "lat"]
    end)
    write3(path, v) = (mkpath(dirname(path)); GREBClimate.jldopen(path, "w") do f
        f["data"] = fill(v, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr); f["dim_names"] = ["lon", "lat", "time"]
    end)
    write_solar(path, v) = (mkpath(dirname(path)); GREBClimate.jldopen(path, "w") do f
        f["data"] = fill(v, GREBClimate.ydim, GREBClimate.nstep_yr); f["dim_names"] = ["lat", "time"]
    end)

    tmpdir = mktempdir()
    try
        write2(joinpath(tmpdir, "static", "global.topography.jld2"), 1.0)
        write2(joinpath(tmpdir, "static", "greb.glaciers.jld2"), 2.0)
        write3(joinpath(tmpdir, "climatology", "ncep.tsurf.1948-2007.clim.jld2"), 3.0)
        write3(joinpath(tmpdir, "climatology", "ncep.zonal_wind.850hpa.clim.jld2"), 4.0)
        write3(joinpath(tmpdir, "climatology", "ncep.meridional_wind.850hpa.clim.jld2"), 5.0)
        write3(joinpath(tmpdir, "climatology", "ncep.atmospheric_humidity.clim.jld2"), 6.0)
        write3(joinpath(tmpdir, "climatology", "ncep.soil_moisture.clim.jld2"), 7.0)
        write3(joinpath(tmpdir, "climatology", "isccp.cloud_cover.clim.jld2"), 8.0)
        write3(joinpath(tmpdir, "climatology", "woce.ocean_mixed_layer_depth.clim.jld2"), 9.0)
        write3(joinpath(tmpdir, "climatology", "Tocean.clim.jld2"), 10.0)
        write3(joinpath(tmpdir, "climatology", "erainterim.omega.vertmean.clim.jld2"), 11.0)
        write3(joinpath(tmpdir, "climatology", "erainterim.omega_std.vertmean.clim.jld2"), 12.0)
        write3(joinpath(tmpdir, "climatology", "erainterim.windspeed.850hpa.clim.jld2"), 13.0)
        write_solar(joinpath(tmpdir, "solar", "solar_radiation.clim.jld2"), 14.0)

        # No flux-corrections file yet: an error that names it, not zeros
        @test_throws "flux_corrections.jld2" load_greb_jld2!(tmpdir; dataset = :ncep)
        @test_throws ArgumentError load_flux_corrections_jld2!(tmpdir, ClimateFields())
        # `corrections = false` loads the rest and leaves the corrections at zero
        without = load_greb_jld2!(tmpdir; dataset = :ncep, corrections = false)
        @test without.loaded && all(==(3.0), without.Tclim)
        @test all(iszero, without.TF_correct) && all(iszero, without.qF_correct) && all(iszero, without.ToF_correct)

        # A dataset name the loader does not know is an error, not NCEP
        @test_throws ArgumentError load_greb_jld2!(tmpdir; dataset = :era5)

        # A file that lacks one of the three tables is an error too
        corrections = joinpath(tmpdir, "climatology", "flux_corrections.jld2")
        GREBClimate.jldopen(corrections, "w") do f
            f["Tsurf_flux_correction"] = fill(15.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
        end
        @test_throws "vapour_flux_correction" load_greb_jld2!(tmpdir; dataset = :ncep)

        GREBClimate.jldopen(corrections, "w") do f
            f["Tsurf_flux_correction"] = fill(15.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
            f["vapour_flux_correction"] = fill(16.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
            f["Tocean_flux_correction"] = fill(17.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
        end

        fields = load_greb_jld2!(tmpdir; dataset = :ncep)
        @test all(==(1.0), fields.z_topo)
        @test all(==(2.0), fields.glacier)
        @test all(==(3.0), fields.Tclim)
        @test all(==(4.0), fields.uclim)
        @test all(==(14.0), fields.sw_solar)
        @test all(==(15.0), fields.TF_correct)
        @test all(==(16.0), fields.qF_correct)
        @test all(==(17.0), fields.ToF_correct)

        # Tasks reading the same file at once all get its content
        path = joinpath(tmpdir, "climatology", "Tocean.clim.jld2")
        reads = [Threads.@spawn read_jld2(path).data for _ in 1:8]
        @test all(t -> all(==(10.0), fetch(t)), reads)

        # A solar table of the wrong shape is reported, with the file's name
        GREBClimate.jldopen(joinpath(tmpdir, "solar", "solar_radiation.clim.jld2"), "w") do f
            f["data"] = zeros(GREBClimate.ydim, 2); f["dim_names"] = ["lat", "time"]
        end
        @test_throws "solar_radiation.clim.jld2" load_greb_jld2!(tmpdir; dataset = :ncep)
    finally
        rm(tmpdir; recursive = true, force = true)
    end

    missing_parent = mktempdir()
    @test_throws ErrorException load_greb_jld2!(joinpath(missing_parent, "nonexistent"))
    rm(missing_parent; recursive = true, force = true)
end

@testset "greb_data_dir resolution order" begin
    tmp_a, tmp_b = mktempdir(), mktempdir()
    saved = get(ENV, "GREB_DATA", nothing)
    try
        # explicit path wins over everything
        ENV["GREB_DATA"] = tmp_b
        @test greb_data_dir(tmp_a) == tmp_a
        # ...and over the environment even with allow_download off
        @test greb_data_dir(tmp_a; allow_download = false) == tmp_a
        # GREB_DATA wins over the repo-local dataset
        @test greb_data_dir() == tmp_b
        delete!(ENV, "GREB_DATA")

        # a non-existent explicit path is an error, not a silent fallback
        @test_throws ErrorException greb_data_dir(joinpath(tmp_a, "nope"))
        # so is a GREB_DATA pointing nowhere
        ENV["GREB_DATA"] = joinpath(tmp_a, "nope")
        @test_throws ErrorException greb_data_dir()
        delete!(ENV, "GREB_DATA")

        # A downloaded dataset is the directory `<load path>/GREB-input-data`.
        # Stand one up so the cache step runs on machines that never downloaded it.
        # DataDeps splits DATADEPS_LOAD_PATH on ':', which cuts a Windows drive
        # letter off, so the load path is given relative to the working directory.
        cd(tmp_b) do
            withenv("DATADEPS_LOAD_PATH" => ".", "DATADEPS_NO_STANDARD_LOAD_PATH" => "true") do
                @test GREBClimate._cached_datadep_path() === nothing
                cache = mkpath(joinpath(pwd(), GREBClimate.DATA_DEP_NAME))
                @test abspath(GREBClimate._cached_datadep_path()) == cache
                # the repo-local dataset wins over the cache; without it the cache is used
                local_dir = normpath(joinpath(@__DIR__, "..", "greb_input_data"))
                @test abspath(greb_data_dir(; allow_download = false)) == (isdir(local_dir) ? local_dir : cache)
            end
        end

        @test greb_data_dir(""; allow_download = false) ==
              greb_data_dir(; allow_download = false)
    finally
        saved === nothing ? delete!(ENV, "GREB_DATA") : (ENV["GREB_DATA"] = saved)
        rm(tmp_a; recursive = true, force = true)
        rm(tmp_b; recursive = true, force = true)
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
             GREBClimate._CC_ANOMALY_FILES, GREBClimate._enso_anomaly_files(:elnino), GREBClimate._FLUX_CORRECTION_KEYS)
    @test all(l -> all(f -> hasfield(ClimateFields, f), keys(l)), lists)

    # The tools keep their own list: it may hold more, never less
    tools = Module()
    Base.include(tools, joinpath(@__DIR__, "..", "tools", "dataset", "fields.jl"))
    @test issubset(names, tools.MODEL_FIELD_NAMES)
    @test issubset(GREBClimate._COMBINED_FILES, tools.COMBINED_FILE_NAMES)
    @test Set(tools.FLUX_CORRECTION_NAMES) == Set(values(GREBClimate._FLUX_CORRECTION_KEYS))
end
