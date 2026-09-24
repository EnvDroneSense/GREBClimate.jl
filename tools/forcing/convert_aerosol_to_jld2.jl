### Convert a GISS stratospheric aerosol NetCDF file to the forcing-series format ###
#
# MAINTAINER TOOL - not part of the package.
#
# Reads one of the zonal-mean files published at
# https://data.giss.nasa.gov/modelforce/strataer/ (Sato-Lacis, CMIP6, CMIP7,
# GloSSAC; all share one layout: 24 latitude bands, monthly steps, `tau` and
# `reff` at 550 nm) and writes a forcing-series file (the layout that
# `_read_forcing_series` in src/io.jl checks) that `load_aerosol_series` reads:
#
#   format_version  Int              2
#   kind            String           "aerosol_aod"
#   calendar        String           "greb_365"
#   time            Vector{Float64}  model decimal years of each month's midpoint
#   lat             Vector{Float64}  band centres, degrees north, ascending
#   values          Matrix{Float64}  stratospheric optical depth at 550 nm (lat x time)
#   units           String           "1"
#   source          String           input file name and its provenance attributes
#   reff            Matrix{Float64}  extra: effective particle radius, micron (lat x time);
#                                    stored for later use, not read by the model yet
#
# Each month is placed at the midpoint of the same month on the model's
# 365-day calendar, so model and record months line up exactly and no date
# conversion happens in the model.
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
#   julia --project=. tools/forcing/convert_aerosol_to_jld2.jl input.nc output.jld2

using NCDatasets
using JLD2
using Dates

const FORMAT_VERSION = 2
const MONTH_DAYS = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)   # the model's 365-day calendar

# Model decimal year of the midpoint of the month containing `t`.
function model_month_midpoint(t)
    m = month(t)
    return year(t) + (sum(MONTH_DAYS[1:m-1]; init=0) + MONTH_DAYS[m] / 2) / 365
end

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
        f["kind"] = "aerosol_aod"
        f["calendar"] = "greb_365"
        f["time"] = model_month_midpoint.(months[keep])
        f["lat"] = lat
        f["values"] = Float64.(tau[:, keep])
        f["units"] = "1"
        f["source"] = source
        f["reff"] = Float64.(reff[:, keep])
    end
    println("wrote $output: $(length(keep)) months, ",
            "$(Dates.format(months[first(keep)], "yyyy-mm")) to $(Dates.format(months[last(keep)], "yyyy-mm")), ",
            "$(length(lat)) latitude bands",
            length(keep) < length(months) ? "; dropped $(length(months) - length(keep)) empty trailing months" : "")
    return output
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 || error("usage: julia --project=. tools/forcing/convert_aerosol_to_jld2.jl input.nc output.jld2")
    convert_aerosol(ARGS[1], ARGS[2])
end
