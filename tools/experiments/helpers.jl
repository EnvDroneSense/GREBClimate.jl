# Shared helpers for experiment scripts: reduce a `greb_model!` result to a few
# numbers and compute the usual diagnostics from maps and series.
#
#   include("helpers.jl"); using .ExperimentTools: summarize, gregory
#
# Pure functions on records, series and maps; none of them runs the model.
# Every mean is area-weighted. The scenario records of a `greb_model!` result are
# anomalies when its `scnr_anomaly` is true; the reducer and the sensitivity
# estimate read that entry.

module ExperimentTools

using GREBClimate

const G = GREBClimate

export annual_means, summarize, in_range_report, gregory, area_mean, zonal_mean, polar_amplification,
       land_ocean, pattern_correlation, sine_fit

const WEIGHT = [Float64(G.dxlat_grid[j]) for _ in 1:G.xdim, j in 1:G.ydim]

"""
    area_mean(map, mask=nothing)

Area-weighted mean of a (`xdim`, `ydim`) map, over the cells where `mask` is
`true` when a mask is given.
"""
function area_mean(map::AbstractMatrix, mask=nothing)
    size(map) == (G.xdim, G.ydim) || throw(DimensionMismatch("expected a ($(G.xdim), $(G.ydim)) map, got $(size(map))"))
    w = mask === nothing ? WEIGHT : WEIGHT .* mask
    sum(w) > 0 || throw(ArgumentError("the mask selects no cell"))
    return sum(Float64.(map) .* w) / sum(w)
end

"The mean over longitude of a map: one value per latitude row."
zonal_mean(map::AbstractMatrix) = vec(sum(Float64.(map); dims=1)) ./ G.xdim

"""
    polar_amplification(map; latitude=60)

`(north, south)`: the mean of `map` poleward of `latitude` degrees over its
global mean, for each hemisphere. A uniform map gives `(1.0, 1.0)`.
"""
function polar_amplification(map::AbstractMatrix; latitude::Real=60)
    lat = G.lat_grid
    rows(mask) = [mask[j] for _ in 1:G.xdim, j in 1:G.ydim]
    g = area_mean(map)
    return (north=area_mean(map, rows(lat .>= latitude)) / g, south=area_mean(map, rows(lat .<= -latitude)) / g)
end

"""
    land_ocean(map, z_topo)

`(land, ocean)`: the area-weighted means of `map` over the land and the ocean
cells, land being where the topography `z_topo` says so.
"""
function land_ocean(map::AbstractMatrix, z_topo::AbstractMatrix)
    land = G.is_land.(z_topo)
    return (land=area_mean(map, land), ocean=area_mean(map, .!land))
end

"Area-weighted correlation of two maps, about their own weighted means."
function pattern_correlation(a::AbstractMatrix, b::AbstractMatrix)
    ma, mb = area_mean(a), area_mean(b)
    da, db = Float64.(a) .- ma, Float64.(b) .- mb
    return sum(WEIGHT .* da .* db) / sqrt(sum(WEIGHT .* da .^ 2) * sum(WEIGHT .* db .^ 2))
end

"""
    sine_fit(series, period)

Least-squares fit of `c + a sin(2 pi t / period) + b cos(2 pi t / period)` to
`series`, `t` in samples from 0. Returns `(amplitude, lag, mean)`: the
amplitude is `hypot(a, b)` and the lag is in samples, positive when the series
peaks later than `sin`.
"""
function sine_fit(series::AbstractVector{<:Real}, period::Real)
    t = 0:(length(series) - 1)
    ω = 2π / period
    design = hcat(ones(length(t)), sin.(ω .* t), cos.(ω .* t))
    c, a, b = design \ Float64.(series)
    # c + A sin(ω (t - lag)) = c + A cos(ω lag) sin(ωt) - A sin(ω lag) cos(ωt)
    lag = atan(-b, a) / ω
    return (amplitude=hypot(a, b), lag=lag, mean=c)
end

"""
    annual_means(records, times, f)

`(years, values)`: for each year with all 12 months among `records`, the mean
over the year of `f(record)`, where `f` is a field name (the area-weighted mean
of that field) or a function of a record that returns a number. `times` holds
the `(year, month)` of each record.
"""
function annual_means(records, times, f)
    length(records) == length(times) || throw(DimensionMismatch("$(length(records)) records, $(length(times)) times"))
    value(rec) = f isa Symbol ? global_mean(getproperty(rec, f)) : Float64(f(rec))
    sums = Dict{Int,Float64}()
    count = Dict{Int,Int}()
    for (rec, t) in zip(records, times)
        sums[t.year] = get(sums, t.year, 0.0) + value(rec)
        count[t.year] = get(count, t.year, 0) + 1
    end
    years = sort!([y for (y, n) in count if n == 12])
    return years, [sums[y] / 12 for y in years]
end

"""
    in_range_report(records, times, anomaly)

`nothing` when every record is finite and inside the allowed range, otherwise
`(year, month, field, value)` of the first one that is not. For absolute
records the range is that of `GREBClimate.RangeCheck()` for `Ts`, `Ta`, `To`
and `q`. For anomalies a value is allowed to be as large as the width of that
range.
"""
function in_range_report(records, times, anomaly::Bool)
    limits = G.RangeCheck().limits
    for (rec, t) in zip(records, times), field in keys(limits)
        lower, upper = limits[field]
        for x in getproperty(rec, field)
            ok = anomaly ? isfinite(x) && abs(x) <= upper - lower : lower <= x <= upper
            ok || return (year=t.year, month=t.month, field=field, value=x)
        end
    end
    return nothing
end

"""
    summarize(result)

Reduce a `greb_model!` result to a few kilobytes: for the control and the
scenario, the annual global means of `Ts`, `precip`, `ice` and the top-of-atmosphere
net (`sw - olr`), with their years; whether the scenario is an anomaly; and
`out_of_range`, the first record outside the range (`nothing` when none; see
[`in_range_report`](@ref)). The control and the scenario are both tested, since a
scenario can leave the range while its control is fine.
"""
function summarize(result)
    series(records, times) = begin
        years, Ts = annual_means(records, times, :Ts)
        (; years, Ts,
           precip=annual_means(records, times, :precip)[2],
           ice=annual_means(records, times, :ice)[2],
           net=annual_means(records, times, rec -> global_mean(rec.sw) - global_mean(rec.olr))[2])
    end
    anomaly = get(result, :scnr_anomaly, false)
    return (ctrl=series(result.ctrl, result.ctrl_time),
            scnr=series(result.scnr, result.scnr_time),
            scnr_anomaly=anomaly,
            out_of_range=(ctrl=in_range_report(result.ctrl, result.ctrl_time, false),
                          scnr=in_range_report(result.scnr, result.scnr_time, anomaly)))
end

"""
    gregory(dTs, dnet)

Gregory estimate from annual global means of a surface-temperature change `dTs`
[K] and a top-of-atmosphere net-flux change `dnet` [W/m2], for example the
`Ts` and `net` series of the scenario in [`summarize`](@ref) of an abrupt-CO2
run. Fits `dnet = forcing + feedback * dTs`.

Returns `forcing` (W/m2), `feedback` (W/m2/K), `sensitivity` (K, the warming at
which the fit reaches zero net flux, `-forcing / feedback`), `r2`, `reached` (the
last `dTs`) and `remaining` (the last `dnet`). The sensitivity can differ from
`reached`: `sw - olr` is not the whole budget of the anomaly.
"""
function gregory(dTs::AbstractVector{<:Real}, dnet::AbstractVector{<:Real})
    length(dTs) == length(dnet) >= 3 || throw(ArgumentError("need at least three annual values of each series"))
    x, y = Float64.(dTs), Float64.(dnet)
    mx, my = sum(x) / length(x), sum(y) / length(y)
    feedback = sum((x .- mx) .* (y .- my)) / sum(abs2, x .- mx)
    forcing = my - feedback * mx
    residual = y .- (forcing .+ feedback .* x)
    r2 = 1 - sum(abs2, residual) / sum(abs2, y .- my)
    return (forcing=forcing, feedback=feedback, sensitivity=-forcing / feedback, r2=r2,
            reached=x[end], remaining=y[end])
end

end # module
