### Convert a GISS stratospheric aerosol NetCDF file to the JLD2 series format ###
#
# MAINTAINER TOOL - not part of the package.
#
# Reads one of the zonal-mean files published at
# https://data.giss.nasa.gov/modelforce/strataer/ (Sato-Lacis, CMIP6, CMIP7,
# GloSSAC; all share one layout: 24 latitude bands, monthly steps, `tau` and
# `reff` at 550 nm) and writes the file `load_aerosol_series` reads:
#
#   format_version  Int              1
#   years           Vector{Float64}  decimal years of the month midpoints
#   lat             Vector{Float64}  band centres, degrees north, ascending
#   aod             Matrix{Float64}  stratospheric optical depth at 550 nm (lat x years)
#   reff            Matrix{Float64}  effective particle radius, micron (lat x years);
#                                    stored for later use, not read by the model yet
#   source          String           input file name and its provenance attributes
#
# Months with no data at all are dropped from the end of the record (the
# Sato-Lacis file pads 2013-2022 with fill values). A gap anywhere else is an
# error, so a truncated or damaged file never becomes a silent zero.
#
# NCDatasets is not a package dependency: add it to your default environment
# once (`julia -e 'using Pkg; Pkg.add("NCDatasets")'`); Julia's environment
# stacking then makes it visible alongside the project.
#
# Usage:
#   julia --project=. tools/convert_aerosol_to_jld2.jl input.nc output.jld2

using NCDatasets
using JLD2
using Dates

const FORMAT_VERSION = 1

decimal_year(t::DateTime) =
    year(t) + (dayofyear(t) - 1 + hour(t) / 24) / daysinyear(t)

function convert_aerosol(input::AbstractString, output::AbstractString)
    lat, months, tau, reff, attrs = NCDataset(input) do ds
        (Float64.(ds["lat"][:]), ds["month"][:], ds["tau"][:, :], ds["reff"][:, :],
         Dict(string(k) => string(v) for (k, v) in ds.attrib))
    end

    # Drop trailing months that are entirely fill values; reject any other gap.
    last_valid = findlast(k -> !all(ismissing, @view tau[:, k]), axes(tau, 2))
    last_valid === nothing && error("$input contains no optical depth data")
    keep = 1:last_valid
    any(ismissing, @view tau[:, keep]) && error("$input has missing optical depths inside the record")
    any(ismissing, @view reff[:, keep]) && error("$input has missing effective radii inside the record")

    issorted(lat) || error("$input: latitudes are not ascending")
    source = join(filter(!isempty, [basename(input), get(attrs, "source", ""),
                                    get(attrs, "history", "")]), " | ")

    jldopen(output, "w") do f
        f["format_version"] = FORMAT_VERSION
        f["years"] = decimal_year.(months[keep])
        f["lat"] = lat
        f["aod"] = Float64.(tau[:, keep])
        f["reff"] = Float64.(reff[:, keep])
        f["source"] = source
    end
    println("wrote $output: $(length(keep)) months, ",
            "$(Dates.format(months[first(keep)], "yyyy-mm")) to $(Dates.format(months[last(keep)], "yyyy-mm")), ",
            "$(length(lat)) latitude bands",
            length(keep) < length(months) ? "; dropped $(length(months) - length(keep)) empty trailing months" : "")
    return output
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 || error("usage: julia --project=. tools/convert_aerosol_to_jld2.jl input.nc output.jld2")
    convert_aerosol(ARGS[1], ARGS[2])
end
