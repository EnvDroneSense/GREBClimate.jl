# The dataset tools (tools/dataset): the CO2 tables, the metadata every file
# carries, the packs and the archives. No raw data needed: small inputs are written here.

# A tool script, loaded into a module of its own (the scripts include each other)
function load_tool(name, file)
    mod = Module(name)
    Core.eval(mod, :(include(path::AbstractString) = Base.include($mod, path)))
    Base.include(mod, joinpath(@__DIR__, "..", "tools", "dataset", file))
    return mod
end
const CONVERTER = load_tool(:DatasetConverter, "convert_greb_to_jld2.jl")
const PACKAGER = load_tool(:DatasetPackager, "package_dataset.jl")

# Run `f()` with what it prints thrown away
function mute(f)
    open(tempname(), "w") do io
        redirect_stdout(f, io)
    end
end

@testset "CO2 tables: rcp60 key, labels, SSP history from the historical record" begin
    with_tempdir() do dir
        hist(y) = 280.0 + 0.5 * (y - 1850)
        write(joinpath(dir, "ipcc.scenario.hist.forcing.CO2.emission.pop.txt"),
            join(("$y $(hist(y)) 2.0 1.5" for y in 1850:2017), "\n") * "\n")
        # RCP tables hold their own history; the SSP tables hold it up to 2008 and the CO2 record from 2009
        rcp(y) = 300.0 + 0.7 * (y - 1950)
        ssp(y) = y < 2009 ? rcp(y) : hist(2009) + 3.0 * (y - 2009) * (y > 2009)
        for key in ("rcp6", "rcp85")
            write(joinpath(dir, "ipcc.scenario.$key.forcing.txt"), join(("$y $(rcp(y))" for y in 1950:2100), "\n") * "\n")
        end
        write(joinpath(dir, "ipcc.scenario.ssp585.forcing.txt"), join(("$y $(ssp(y))" for y in 1950:2100), "\n") * "\n")

        out = joinpath(dir, "out")
        mute(() -> CONVERTER.convert_scenario_texts(dir, joinpath(out, "scenario")))
        file = joinpath(out, "scenario", "ipcc_scenarios.jld2")
        scenarios = jldopen(f -> f["scenarios"], file)
        labels = jldopen(f -> f["labels"], file)

        @test sort(collect(keys(scenarios))) == ["hist", "rcp60", "rcp85", "ssp585"]   # "rcp6" is "rcp60"
        @test labels == Dict("hist" => "co2", "rcp60" => "co2_equivalent", "rcp85" => "co2_equivalent", "ssp585" => "co2")
        # RCP tables are as in the raw files
        @test all(y -> scenarios["rcp85"][y] ≈ rcp(y), 1950:2100)
        @test scenarios["rcp60"] == scenarios["rcp85"]
        # The SSP table: the CO2 record before 2009, its own values from 2009 on
        @test all(y -> scenarios["ssp585"][y] ≈ hist(y), 1950:2008)
        @test all(y -> scenarios["ssp585"][y] ≈ ssp(y), 2009:2100)
        @test length(scenarios["ssp585"]) == 151
        # The file reads back through the package
        @test load_co2_scenario(out, :ssp585)[1950] ≈ Float32(hist(1950))
        @test co2_scenario_label(out, :rcp60) === :co2_equivalent

        # An SSP table that does not meet the CO2 record in 2009 is not what the splice assumes
        write(joinpath(dir, "ipcc.scenario.ssp585.forcing.txt"), join(("$y $(rcp(y))" for y in 1950:2100), "\n") * "\n")
        @test_throws "do not meet" mute(() -> CONVERTER.convert_scenario_texts(dir, joinpath(dir, "out2")))
    end
end

@testset "every file carries its metadata" begin
    with_tempdir() do dir
        info = CONVERTER.field_info("global.topography")
        path = joinpath(dir, "static", "global.topography.jld2")
        CONVERTER.write_jld2_field(path, zeros(Float32, 96, 48), ["lon", "lat"], info)
        field = read_field(path)
        @test field.info["units"] == "m" && field.info["dataset_version"] == GREBClimate.DATASET_VERSION
        @test all(k -> haskey(field.info, k), ("description", "source", "period"))
        # The axes: longitude from 0, latitude from -88.125, both in steps of 3.75
        @test field.coords[1][[1, 96]] == [0.0, 356.25] && field.coords[2][[1, 48]] == [-88.125, 88.125]
        @test PACKAGER.metadata_problems(path, "global.topography") == String[]

        # A file without a key, with another version, or with axis values of the wrong length is reported
        jldopen(path, "w") do file
            file["data"] = zeros(Float32, 96, 48)
            file["dim_names"] = ["lon", "lat"]
            file["coords"] = Dict(1 => [0.0, 3.75])
            file["dataset_version"] = 1
        end
        problems = PACKAGER.metadata_problems(path, "global.topography")
        @test any(contains("lacks \"units\""), problems)
        @test any(contains("dataset_version 1"), problems)
        @test any(contains("axis 1 has 2 values for 96"), problems)

        # Every field the model reads has metadata, and a name without any is an error, not a blank
        @test all(name -> !isempty(CONVERTER.field_info(name).source), CONVERTER.MODEL_FIELD_NAMES)
        @test_throws ErrorException CONVERTER.field_info("unknown.field")
    end
end

@testset "packager: the tree is validated, the archive is reproducible" begin
    with_tempdir() do dir
        # An empty tree lacks every file; a stray file is refused
        @test_throws "missing" PACKAGER.validate_dataset(dir)
        @test_throws ErrorException PACKAGER.validate_dataset(joinpath(dir, "nope"))

        if Sys.which("tar") === nothing || Sys.which("xz") === nothing
            @test_skip "tar and xz are needed to build archives"
        else
            tree = joinpath(dir, "tree")
            mkpath(joinpath(tree, "static"))
            write(joinpath(tree, "static", "a.jld2"), rand(UInt8, 5000))
            write(joinpath(tree, "static", "b.jld2"), rand(UInt8, 5000))
            files = PACKAGER.dataset_files(tree)
            @test files == ["static/a.jld2", "static/b.jld2"]

            one = PACKAGER.build_archive(tree, files, joinpath(dir, "one.tar.xz"))
            # Same contents, later timestamps: the same bytes
            sleep(1.1)
            foreach(f -> touch(joinpath(tree, f)), files)
            two = PACKAGER.build_archive(tree, files, joinpath(dir, "two.tar.xz"))
            @test read(one) == read(two)

            # The archive holds exactly the listed files
            unpacked = mkpath(joinpath(dir, "unpacked"))
            run(pipeline(`xz -dc $one`, `tar -xf - -C $unpacked`))
            @test PACKAGER.dataset_files(unpacked) == files
            @test read(joinpath(unpacked, "static", "a.jld2")) == read(joinpath(tree, "static", "a.jld2"))
        end
    end
    @test PACKAGER.archive_name(:core) == "greb_input_data-v$(GREBClimate.DATASET_VERSION)-core.tar.xz"
end

@testset "the tools' file list, packs and version agree with the loaders" begin
    names = GREBClimate.dataset_field_files()

    # The tools keep their own list: it may hold more, never less
    tools = Module()
    Base.include(tools, joinpath(@__DIR__, "..", "tools", "dataset", "fields.jl"))
    @test issubset(names, tools.MODEL_FIELD_NAMES)
    @test issubset(GREBClimate._COMBINED_FILES, tools.COMBINED_FILE_NAMES)
    @test Set(tools.FLUX_CORRECTION_NAMES) == Set(values(GREBClimate._FLUX_CORRECTION_KEYS))

    # The tools write the version the loaders read
    @test tools.DATASET_VERSION == GREBClimate.DATASET_VERSION
    # The tools put a file in the pack the loaders name for it, and nothing else in an optional pack
    @test Tuple(tools.PACK_NAMES) == keys(GREBClimate.DATA_PACKS)
    optional = GREBClimate._OPTIONAL_PACK_FILES
    for pack in GREBClimate._OPTIONAL_PACKS, file in getfield(optional, pack)

        @test tools.pack_of(file) === pack
    end
    in_core = [f for f in ("static/global.topography.jld2", "solar/solar_radiation.clim.jld2",
        "climatology/ncep.tsurf.1948-2007.clim.jld2", "climatology/flux_corrections.jld2",
        "climatology/erainterim.omega.vertmean.clim.jld2", "scenario/ipcc_scenarios.jld2")]
    @test all(f -> tools.pack_of(f) === :core, in_core)
    # Every single-field file belongs to exactly one pack, and the tools have metadata for each
    @test all(name -> tools.field_info(name).units isa String, tools.MODEL_FIELD_NAMES)
end
