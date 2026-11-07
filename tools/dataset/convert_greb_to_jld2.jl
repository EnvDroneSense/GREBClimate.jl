# Convert the raw GREB files (layout in DATA_README.md) into the `.jld2` dataset
# read by `load_climatology` and `GREBClimate.read_field`. MAINTAINER TOOL: users
# get the prepared packs (see README.md).
#
# Each field becomes one file with the keys "data", "dim_names", "coords" (the
# values of each axis), "description", "units", "source", "period" and
# "dataset_version"; a GrADS `.bin` also keeps the text of its `.ctl` under "ctl".
# The flux corrections, the solar scenario families and the CO2 text files are
# combined into one file each (see below). Any failure stops the script.
#
#   julia --project=. tools/dataset/convert_greb_to_jld2.jl [input_dir] [output_dir] [--all]
#   # defaults: Data, greb_input_data
#   # --all also converts fields the model never reads (the packager refuses them)

using JLD2

include(joinpath(@__DIR__, "fields.jl"))   # file list, packs, metadata text

const LON, LAT, TIME = 96, 48, 730

# the layout of a raw file, by its size in bytes
const LAYOUTS = Dict(
    LON * LAT * 4 => ((LON, LAT), ["lon", "lat"]),
    LON * LAT * TIME * 4 => ((LON, LAT, TIME), ["lon", "lat", "time"]),
    LAT * TIME * 4 => ((LAT, TIME), ["lat", "time"])
)

function read_bin(path::AbstractString, dims::Tuple)
    data = Vector{Float32}(undef, prod(dims))
    read!(path, data)
    return reshape(data, dims)
end

function read_ctl_text(bin_path::AbstractString)
    ctl_path = splitext(bin_path)[1] * ".ctl"
    return isfile(ctl_path) ? read(ctl_path, String) : nothing
end

# Axis values of the axes named in `dim_names`, keyed by axis number
axis_coords(dim_names) = Dict{Int,Vector{Float64}}(
    i => GRID_COORDS[Symbol(n)] for (i, n) in enumerate(dim_names) if haskey(GRID_COORDS, Symbol(n)))

# The metadata of `info` (description, units, source, and period when it has one)
function write_info(file, info)
    for key in (:description, :units, :source, :period)
        haskey(info, key) && (file[String(key)] = info[key])
    end
    file["dataset_version"] = DATASET_VERSION
end

"Write a field file: the array, its axes, the metadata in `info`. `coords` adds or replaces axis values."
function write_jld2_field(path::AbstractString, data, dim_names::Vector{String}, info;
        coords::AbstractDict=Dict{Int,Vector{Float64}}(), ctl::Union{Nothing,String}=nothing)
    mkpath(dirname(path))
    jldopen(path, "w") do file
        file["data"] = data
        file["dim_names"] = dim_names
        file["coords"] = merge(axis_coords(dim_names), coords)
        write_info(file, info)
        ctl === nothing || (file["ctl"] = ctl)
    end
end

# ---- the single fields ----

function convert_fields(input_dir, output_dir; all::Bool)
    bins = filter(p -> endswith(p, ".bin") && splitext(basename(p))[1] ∉ FLUX_CORRECTION_NAMES, readdir(input_dir; join=true))
    name_of(p) = splitext(basename(p))[1]
    skipped = filter(p -> name_of(p) ∉ MODEL_FIELD_NAMES, bins)
    all || filter!(p -> name_of(p) ∈ MODEL_FIELD_NAMES, bins)
    all || isempty(skipped) || println("  skipped $(length(skipped)) raw files the model never reads (--all converts them)")

    absent = setdiff(MODEL_FIELD_NAMES, name_of.(bins))
    isempty(absent) || error("$input_dir lacks $(length(absent)) field(s) the model reads: ", join(sort(collect(absent)), ", "))

    folders = String[]
    for path in sort(bins)
        name = name_of(path)
        dims, dim_names = LAYOUTS[stat(path).size]
        folder = occursin(r"topography|glacier", name) ? "static" : occursin("solar_radiation", name) ? "solar" : "climatology"
        write_jld2_field(joinpath(output_dir, folder, "$name.jld2"), read_bin(path, dims), dim_names, field_info(name);
            ctl=read_ctl_text(path))
        push!(folders, folder)
    end
    println("  fields: $(length(bins)) raw .bin files -> ", join(("$(count(==(f), folders)) in $f/" for f in unique(folders)), ", "))
end

"The three flux corrections in one file, keyed by field name."
function convert_flux_corrections(input_dir, output_dir)
    path = joinpath(output_dir, "climatology", "flux_corrections.jld2")
    mkpath(dirname(path))
    jldopen(path, "w") do file
        for name in FLUX_CORRECTION_NAMES
            file[name] = read_bin(joinpath(input_dir, "$name.bin"), (LON, LAT, TIME))
        end
        dim_names = ["lon", "lat", "time"]
        file["dim_names"] = dim_names
        file["coords"] = axis_coords(dim_names)
        file["description"] = Dict(String(k) => v.description for (k, v) in pairs(FLUX_CORRECTION_INFO))
        file["units"] = Dict(String(k) => v.units for (k, v) in pairs(FLUX_CORRECTION_INFO))
        write_info(file, (source=SOURCE_GREB, period="mean year"))
    end
    println("  flux corrections: ", join(FLUX_CORRECTION_NAMES, ", "), " -> climatology/flux_corrections.jld2")
end

# ---- the solar scenarios: the paleo file, and the eccentricity and obliquity families ----

"Stack the files `greb.solar.<prefix>.<N>.bin` of `dir`, sorted by the number N, along a new first axis called `index_name`."
function convert_solar_family(dir, prefix, out_path, index_name, description)
    found = [(parse(Float64, m[1]), joinpath(dir, f)) for f in readdir(dir)
             for m in (match(Regex("^greb\\.solar\\.$prefix\\.(\\d+)\\.bin\$"), f),) if m !== nothing]
    isempty(found) && error("no files greb.solar.$prefix.<N>.bin in $dir")
    sort!(found; by=first)
    stack = permutedims(cat((read_bin(p, (LAT, TIME)) for (_, p) in found)...; dims=3), (3, 1, 2))
    write_jld2_field(out_path, stack, [index_name, "lat", "time"],
        (description=description, units="W m-2", source=SOURCE_GREB, period="mean year");
        coords=Dict(1 => first.(found), 2 => GRID_COORDS.lat, 3 => GRID_COORDS.time))
    return length(found)
end

function convert_solar_scenarios(input_dir, output_dir)
    dir = joinpath(input_dir, "solar_forcing_scenarios")
    isdir(dir) || error("$dir not found (the solar pack needs it)")
    paleo = joinpath(dir, "greb.solar.231K_hybers.corrected.bin")
    write_jld2_field(joinpath(output_dir, "solar_paleo.jld2"), read_bin(paleo, (LAT, TIME)), ["lat", "time"],
        (description="Top-of-atmosphere solar radiation of the 231 kyr paleo experiment", units="W m-2",
            source=SOURCE_GREB, period="mean year"); ctl=read_ctl_text(paleo))
    counts = map((("eccentricity", "ecc_index"), ("obliquity", "obl_index"))) do (prefix, index_name)
        convert_solar_family(dir, prefix, joinpath(output_dir, "solar_$prefix.jld2"), index_name,
            "Top-of-atmosphere solar radiation for each orbital $prefix scenario")
    end
    println("  solar scenarios: the paleo table, $(counts[1]) eccentricity and $(counts[2]) obliquity tables -> solar_scenarios/")
end

# ---- the CO2 text files (RCP, SSP, historical), combined into one file ----

# The Fortran file names call the RCP6.0 table "rcp6"; the dataset calls it "rcp60".
const TABLE_KEY_RENAMES = Dict("rcp6" => "rcp60")

# The SSP files hold the RCP (CO2-equivalent) history up to 2008 and CO2 only from
# this year on. The dataset replaces their earlier years with the historical CO2
# record, so each SSP table is CO2 only throughout.
const SSP_FIRST_CO2_YEAR = 2009

"What a table holds: `co2` for the SSP and historical tables, `co2_equivalent` for the RCP tables."
co2_label(key::AbstractString) = startswith(key, "rcp") ? "co2_equivalent" : "co2"

"Read a text file of `year value ...` lines into a `year => column` table (`column` is the column of the value)."
function parse_table(path::AbstractString, column::Int=2)
    cols = [split(strip(line)) for line in eachline(path) if !isempty(strip(line))]
    return Dict(parse(Int, c[1]) => parse(Float64, c[column]) for c in cols)
end

"""
Replace the years before `SSP_FIRST_CO2_YEAR` of every SSP table with the historical
CO2 record `hist`, and return the set of years replaced. An SSP table must equal `hist`
in its first CO2-only year, or the tables are not what this assumes.
"""
function splice_co2_history!(scenarios::Dict{String,Dict{Int,Float64}})
    hist = get(scenarios, "hist") do
        error("the historical CO2 table (ipcc.scenario.hist.forcing*.txt) is missing")
    end
    replaced = Set{Int}()
    for (key, table) in scenarios
        startswith(key, "ssp") || continue
        isapprox(table[SSP_FIRST_CO2_YEAR], hist[SSP_FIRST_CO2_YEAR]; atol=0.01) ||
            error("$key in $SSP_FIRST_CO2_YEAR is $(table[SSP_FIRST_CO2_YEAR]) but the historical CO2 is " *
                  "$(hist[SSP_FIRST_CO2_YEAR]): the tables do not meet")
        for year in filter(<(SSP_FIRST_CO2_YEAR), collect(keys(table)))
            table[year] = get(hist, year) do
                error("the historical CO2 table has no value for $year (needed by $key)")
            end
            push!(replaced, year)
        end
    end
    return replaced
end

"""
Combine every `ipcc.scenario.<key>.forcing*.txt` of `input_dir` into
`ipcc_scenarios.jld2` (RCP6.0 as `rcp60`, with a label per table), and the emissions
and population columns of the historical file into `historical_emissions_population.jld2`.
"""
function convert_scenario_texts(input_dir, output_dir)
    pattern = r"^ipcc\.scenario\.(.+)\.forcing.*\.txt$"
    paths = [(String(get(TABLE_KEY_RENAMES, m[1], m[1])), joinpath(input_dir, f))
             for f in readdir(input_dir) for m in (match(pattern, f),) if m !== nothing]
    scenarios = Dict(key => parse_table(path) for (key, path) in paths)
    replaced = splice_co2_history!(scenarios)

    mkpath(output_dir)
    jldopen(joinpath(output_dir, "ipcc_scenarios.jld2"), "w") do file
        file["scenarios"] = scenarios
        file["labels"] = Dict(key => co2_label(key) for key in keys(scenarios))
        file["periods"] = Dict(key => "$(minimum(keys(t)))-$(maximum(keys(t)))" for (key, t) in scenarios)
        write_info(file, (description="CO2 per year of each scenario. labels: co2 is CO2 only, co2_equivalent is an " *
                                      "equivalent concentration that includes the other gases and aerosol. " *
                                      "The SSP tables hold the historical CO2 record before $SSP_FIRST_CO2_YEAR.",
            units="ppm", source=SOURCE_IPCC))
    end

    # columns 3-4 of the historical file: CO2 emissions in Gt CO2/yr, world population in billions
    hist_path = last(paths[findfirst(p -> first(p) == "hist", paths)])
    emissions = parse_table(hist_path, 3)
    population = parse_table(hist_path, 4)
    data = Dict(year => (co2_emissions_gt_co2_yr=emissions[year], population_billions=population[year])
    for year in sort(collect(keys(emissions))))
    jldopen(joinpath(output_dir, "historical_emissions_population.jld2"), "w") do file
        file["data"] = data
        write_info(file, (description="Historical CO2 emissions and world population per year",
            units=Dict("co2_emissions_gt_co2_yr" => "Gt CO2 yr-1", "population_billions" => "1e9"),
            source=SOURCE_IPCC, period="$(minimum(keys(data)))-$(maximum(keys(data)))"))
    end
    ssp = sort!([k for k in keys(scenarios) if startswith(k, "ssp")])
    println("  CO2 tables: ", join(sort!(collect(keys(scenarios))), ", "), " -> scenario/ipcc_scenarios.jld2")
    println("    rcp6 is stored as rcp60; $(join(ssp, ", ")) take the historical CO2 for ",
        "$(minimum(replaced))-$(maximum(replaced)) (labels: co2 for ssp* and hist, co2_equivalent for rcp*)")
    println("  historical emissions and population, $(minimum(keys(data)))-$(maximum(keys(data))) -> scenario/historical_emissions_population.jld2")
end

function main(input_dir, output_dir; all::Bool=false)
    isdir(input_dir) || error("""
        Raw input directory not found: $input_dir (layout: DATA_README.md). Pass it explicitly:
            julia --project=. tools/dataset/convert_greb_to_jld2.jl <input_dir> [output_dir]""")
    println("Converting $input_dir -> $output_dir")
    convert_fields(input_dir, output_dir; all)
    convert_flux_corrections(input_dir, output_dir)
    convert_solar_scenarios(input_dir, joinpath(output_dir, "solar_scenarios"))
    convert_scenario_texts(input_dir, joinpath(output_dir, "scenario"))
    files = [f for (root, _, names) in walkdir(output_dir) for f in joinpath.(root, names) if endswith(f, ".jld2")]
    println("Done: $(length(files)) .jld2 files, $(round(sum(filesize, files) / 1048576; digits=1)) MB in $output_dir.")
    println("Next: julia --project=. tools/dataset/package_dataset.jl $output_dir (checks the tree, builds the packs).")
end

if abspath(PROGRAM_FILE) == @__FILE__
    positional = filter(a -> !startswith(a, "--"), ARGS)
    main(length(positional) >= 1 ? positional[1] : joinpath(@__DIR__, "..", "..", "Data"),
        length(positional) >= 2 ? positional[2] : joinpath(@__DIR__, "..", "..", "greb_input_data");
        all="--all" in ARGS)
end
