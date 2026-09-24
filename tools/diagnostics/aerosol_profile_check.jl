### Aerosol kernel against the published records, per latitude band ###
#
# MAINTAINER TOOL - not part of the package. Needs the converted aerosol
# records in Data/aerosol/ (tools/forcing/convert_aerosol_to_jld2.jl); no model
# run and no dataset.
#
# For each reference eruption, runs `aerosol_optical_depth!` with the
# eruption's class, scaled so its global-mean peak equals the record's, and
# compares it with the record on the model's 48 latitude rows:
#
#   band ratio    each band's peak optical depth over the largest band's
#   peak month    months from the eruption date to each band's maximum
#   1/e months    months from each band's maximum until it falls to 1/e of it
#                 (both left blank for bands below 5 percent of the largest)
#   NH/SH         northern over southern area-mean optical depth, summed over the window
#   shape error   RMS difference of the global-mean curves, each scaled to a peak of 1
#
# The record's background (its mean over the 12 months before the eruption, or
# before `background` where an earlier eruption is still decaying) is subtracted
# per latitude row, as in tools/validation/validate_volcanic.jl.
#
# Katmai is also compared with Stothers (1996) Table 6, the reference for the
# extratropical classes; see `STOTHERS_KATMAI`.
# Plan: aerosol-parameter-plan (vault), tasks 1 and 1b.
#
# Usage:
#   julia --project=. tools/diagnostics/aerosol_profile_check.jl

using GREBClimate
using Printf

const AEROSOL_DIR = normpath(joinpath(@__DIR__, "..", "..", "Data", "aerosol"))
const SHAPE_MONTHS = 36          # months used for the shape error
const LAT = Float64.(GREBClimate.lat_grid)
const W = cosd.(LAT)

# Latitude bands per hemisphere, degrees from the equator.
const BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 90))
const BAND_NAMES = [string(lo, "-", hi, h > 0 ? "N" : "S") for h in (1, -1) for (lo, hi) in BANDS]
const ROWS = [findall(φ -> lo <= h * φ < hi, LAT) for h in (1, -1) for (lo, hi) in BANDS]
const NH, SH = findall(>(0), LAT), findall(<(0), LAT)

# The class a reference eruption should use, falling back while the class is
# not yet implemented (one-hemisphere classes: plan task 4).
class_or_tropical(cls) = cls in GREBClimate.AEROSOL_CLASSES ? cls : :tropical

# `window`: months compared, cut short where the next eruption begins.
# `background`: end of the 12 months used as background (default: the eruption date).
# `note`: printed with the result.
const ERUPTIONS = (
    (name="Pinatubo", date=1991.45, class=:tropical, records=("Sato-Lacis", "GloSSAC.v2.24"), window=48),
    (name="El Chichon", date=1982.26, class=:tropical_nh, records=("Sato-Lacis", "GloSSAC.v2.24"), window=48),
    (name="Agung", date=1963.21, class=:tropical_sh, records=("Sato-Lacis",), window=48),
    (name="Katmai", date=1912.43, class=:nh_extratropical, records=("Sato-Lacis",), window=48,
     note="Sato-Lacis before 1960 spreads one station per band: Simla (31N) for 0-30N, Pavlovsk (56N) for 30-90N (Sato et al. 1993)"),
    (name="Kasatochi", date=2008.60, class=:nh_extratropical, records=("GloSSAC.v2.24",), window=10,
     note="report only: about 30 times smaller than Katmai and injected lower; window ends before Sarychev"),
    (name="Sarychev", date=2009.45, class=:nh_extratropical, records=("GloSSAC.v2.24",), window=16,
     background=2008.60, note="report only: as Kasatochi; background taken before Kasatochi"),
)

# Stothers (1996) Table 6: pyrheliometric optical depth x 1000 after Katmai
# (58 N, 6 June 1912), June 1912 to September 1914. Blank entries in the table
# are zero; there are none south of 30 N. The month columns are fixed by the
# counts: 7 values in 1912 (June to December), 12 in 1913, 9 in 1914.
const STOTHERS_KATMAI = (
    date = 1912.43,
    start = (1912, 6),
    bands = ((30, 45), (45, 90)),
    values = (
        [29, 86, 112, 95, 84, 66, 63,
         70, 83, 75, 74, 39, 55, 30, 31, 19, 33, 23, 26,
         43, 28, 21, 3, 20, 9, 0, 21, 0],
        [40, 237, 223, 228, 118, 115, 76,
         69, 52, 65, 71, 65, 44, 58, 49, 18, 23, 50, 37,
         60, 33, 25, 14, 31, 8, 46, 27, 12],
    ),
    efold_years = (0.8, 0.1),    # Stothers' own least-squares fit
)

area_mean(v, rows) = sum(v[rows] .* W[rows]) / sum(W[rows])

# Record times and latitude-by-month optical depth after `date`, background removed.
function record_window(series, date, window, background)
    k0 = findfirst(>=(date), series.years)
    before = findall(t -> background - 1 <= t < background, series.years)
    (k0 === nothing || isempty(before) || k0 + window - 1 > length(series.years)) && return nothing
    bg = sum(series.aod[:, before]; dims=2) ./ length(before)
    return series.years[k0:k0 + window - 1], Float64.(series.aod[:, k0:k0 + window - 1] .- bg)
end

function model_window(times, date, cls, peak)
    sc = AerosolScenario(eruptions=[Eruption(date, cls; peak_aod=peak)])
    buf = zeros(Float32, ydim)
    return reduce(hcat, (GREBClimate.aerosol_optical_depth!(buf, sc, t); Float64.(buf)) for t in times)
end

# Months from a series' maximum until it falls to 1/e of it; `missing` if it never does.
function efold_months(s, months)
    i = argmax(s)
    j = findfirst(<=(s[i] / exp(1)), s[i:end])
    return j === nothing ? missing : months[i + j - 1] - months[i]
end

# e-folding time in months from a least-squares fit of log(s) after the maximum,
# the method Stothers (1996) used; non-positive values are skipped.
function efold_fit(s, months)
    k = [j for j in argmax(s):length(s) if s[j] > 0]
    x, y = months[k], log.(s[k])
    slope = sum((x .- sum(x) / length(x)) .* (y .- sum(y) / length(y))) / sum(abs2, x .- sum(x) / length(x))
    return -1 / slope
end

# `months` holds each column's time since the eruption, in months.
function metrics(aod, months)
    series = [[area_mean(aod[:, k], r) for k in axes(aod, 2)] for r in ROWS]
    ratio = maximum.(series) ./ maximum(maximum.(series))
    visible = ratio .>= 0.05
    peak_month = [v ? months[argmax(s)] : nothing for (s, v) in zip(series, visible)]
    efold = [v ? efold_months(s, months) : nothing for (s, v) in zip(series, visible)]
    gm = [area_mean(aod[:, k], eachindex(LAT)) for k in axes(aod, 2)]
    nh, sh = (sum(area_mean(aod[:, k], h) for k in axes(aod, 2)) for h in (NH, SH))
    nhsh = nh > 0 && sh >= 0 ? nh / sh : NaN    # undefined when a hemisphere is below background
    return (ratio=ratio, peak_month=peak_month, efold=efold, gm=gm,
            gm_peak=months[argmax(gm)], nhsh=nhsh)
end

cell(::Nothing) = "    -"
cell(::Missing) = @sprintf("%5s", ">win")
cell(x::Real) = @sprintf("%5.1f", x)

function report(e, cls, recname, series)
    win = record_window(series, e.date, e.window, get(e, :background, e.date))
    win === nothing && return println("\n$(e.name) / $recname: record does not cover the window")
    times, rec = win
    peak = maximum(area_mean(rec[:, k], eachindex(LAT)) for k in axes(rec, 2))
    months = (times .- e.date) .* 12
    r = metrics(rec, months)
    m = metrics(model_window(times, e.date, cls, peak), months)
    n = min(SHAPE_MONTHS, e.window)
    shape = sqrt(sum(abs2, r.gm[1:n] ./ maximum(r.gm) .- m.gm[1:n] ./ maximum(m.gm)) / n)
    @printf("\n%s (%.2f), class %s, record %s, %d months: global-mean peak %.4f; NH/SH record %.2f, model %.2f; shape error %.3f\n",
            e.name, e.date, cls, recname, e.window, peak, r.nhsh, m.nhsh, shape)
    haskey(e, :note) && println("  ($(e.note))")
    println("  band       ratio rec   mod   peak month rec   mod   1/e months rec   mod")
    for (i, b) in enumerate(BAND_NAMES)
        @printf("  %-8s       %5.2f %5.2f              %s %s              %s %s\n", b, r.ratio[i], m.ratio[i],
                cell(r.peak_month[i]), cell(m.peak_month[i]), cell(r.efold[i]), cell(m.efold[i]))
    end
    @printf("  %-8s                                %s %s\n", "global", cell(r.gm_peak), cell(m.gm_peak))
end

# Katmai against Stothers (1996) Table 6, on its two bands. Ratios and times
# only: pyrheliometric values are not on the same scale as the model's
# (Stothers: visual = 1.6 x pyrheliometric).
function report_stothers(cls)
    s = STOTHERS_KATMAI
    n = length(s.values[1])
    times = [s.start[1] + (s.start[2] - 1 + k + 0.5) / 12 for k in 0:n - 1]
    months = (times .- s.date) .* 12
    aod = model_window(times, s.date, cls, 1.0)
    rows(lo, hi) = findall(φ -> lo <= φ < hi, LAT)
    rec = [Float64.(v) for v in s.values]
    mdl = [[area_mean(aod[:, k], rows(lo, hi)) for k in 1:n] for (lo, hi) in s.bands]
    tropics = maximum(area_mean(aod[:, k], rows(0, 30)) for k in 1:n) / maximum(maximum.(mdl))
    @printf("\nKatmai (%.2f), class %s, against Stothers (1996) Table 6, %d months\n", s.date, cls, n)
    println("  band       ratio rec   mod   peak month rec   mod   1/e months rec   mod   fitted e-fold rec   mod")
    top_r, top_m = maximum(maximum.(rec)), maximum(maximum.(mdl))
    for (i, (lo, hi)) in enumerate(s.bands)
        @printf("  %-8s       %5.2f %5.2f              %s %s              %s %s                 %s %s\n",
                "$lo-$(hi)N", maximum(rec[i]) / top_r, maximum(mdl[i]) / top_m,
                cell(months[argmax(rec[i])]), cell(months[argmax(mdl[i])]),
                cell(efold_months(rec[i], months)), cell(efold_months(mdl[i], months)),
                cell(efold_fit(rec[i], months)), cell(efold_fit(mdl[i], months)))
    end
    @printf("  0-30N    peak over the largest band: Stothers 0 (no aerosol at Helwan 30N, Mexico City 19N, Arequipa 16S), model %.2f\n", tropics)
    @printf("  Stothers' published e-folding time: %.1f +/- %.1f months\n", 12 * s.efold_years[1], 12 * s.efold_years[2])
end

function main()
    series = Dict(n => load_aerosol_series(joinpath(AEROSOL_DIR, "aerosol_$n.jld2"))
                  for n in ("Sato-Lacis", "GloSSAC.v2.24"))
    for e in ERUPTIONS
        cls = class_or_tropical(e.class)
        cls == e.class || println("\n(note: $(e.name) uses :$cls; :$(e.class) is not implemented yet)")
        for rec in e.records
            report(e, cls, rec, series[rec])
        end
        e.name == "Katmai" && report_stothers(cls)
    end
end

main()
