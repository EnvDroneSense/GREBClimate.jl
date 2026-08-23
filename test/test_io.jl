# JLD2 loading, dataset resolution, and converter/archive consistency.

@testset "read_jld2 rejects non-JLD2 input" begin
    tmp = tempname() * ".jld2"
    write(tmp, "not a jld2 file")
    @test_throws Exception read_jld2(tmp)
    rm(tmp; force = true)
end

@testset "load_greb_jld2!/load_flux_corrections_jld2! file-exists branches" begin
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

        # "files missing" branch: no flux-correction files present yet ->
        # load_flux_corrections_jld2! should warn and zero-fill, not error.
        fields_nocorr = load_greb_jld2!(tmpdir; dataset = :ncep)
        @test all(==(0.0), fields_nocorr.TF_correct)
        @test all(==(0.0), fields_nocorr.qF_correct)
        @test all(==(0.0), fields_nocorr.ToF_correct)
        @test all(==(3.0), fields_nocorr.Tclim)  # loader itself still worked

        # "files present" branch: add the combined flux-correction file and reload.
        mkpath(joinpath(tmpdir, "climatology"))
        GREBClimate.jldopen(joinpath(tmpdir, "climatology", "flux_corrections.jld2"), "w") do f
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

        cached = GREBClimate._cached_datadep_path()
        if cached === nothing
            @test_skip "no DataDeps cache on this machine"
        else
            @test isdir(cached)
            @test greb_data_dir(; allow_download = false) !== nothing
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

@testset "converter allowlist matches what src/io.jl loads" begin
    repo = normpath(joinpath(@__DIR__, ".."))
    conv = read(joinpath(repo, "tools", "convert_greb_to_jld2.jl"), String)
    io_src = read(joinpath(repo, "src", "io.jl"), String)

    # --- the allowlist, as literals inside the MODEL_FIELD_NAMES block ---
    m = match(r"const MODEL_FIELD_NAMES = Set\{String\}\(\[(.*?)
\]\)"s, conv)
    @test m !== nothing
    # Cut at the ENSO comprehension: its "zonal.wind"/"meridional.wind"
    # tokens are field-name *fragments*, not file names, and it is expanded
    # explicitly below.
    body = m.captures[1]
    cut = findfirst("(\"erainterim.", body)
    cut === nothing || (body = body[1:first(cut)-1])
    allowed = Set{String}()
    for lit in eachmatch(r"\"([^\"]+)\"", body)
        s = lit.captures[1]
        if occursin('$', s)
            continue          # the ENSO comprehension template, expanded below
        elseif occursin('.', s) && !occursin(' ', s)
            push!(allowed, s)
        end
    end
    # expand the ENSO comprehension the same way the converter does
    for f in ("tsurf", "zonal.wind", "meridional.wind", "windspeed", "omega"),
        s in ("elnino", "lanina")
        push!(allowed, "erainterim.$f.$s.forcing")
    end
    @test length(allowed) == 33

    # --- what io.jl actually loads, with $suffix expanded ---
    loaded = Set{String}()
    for m2 in eachmatch(r"\"([A-Za-z0-9_.\$-]+)\.jld2\"", io_src)
        name = m2.captures[1]
        # combined multi-field files are not per-field entries in the allowlist
        name in ("flux_corrections", "ipcc_scenarios", "solar_paleo",
                 "solar_eccentricity", "solar_obliquity") && continue
        if occursin("\$suffix", name)
            for s in ("elnino", "lanina")
                push!(loaded, replace(name, "\$suffix" => s))
            end
        else
            push!(loaded, name)
        end
    end

    # Every field io.jl loads must be produced by the converter, and the
    # converter must not carry entries nothing loads.
    @test isempty(setdiff(loaded, allowed))
    @test isempty(setdiff(allowed, loaded))
end
