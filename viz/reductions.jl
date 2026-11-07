# Records in, numbers out. Records are NamedTuples of lon x lat matrices; the
# grid size is read from the data, so nothing here imports the model.

"Cell-centre latitudes [°N] of an `ny`-row grid (the model's `lat_grid`)."
lats(ny::Integer) = [180 / ny * (k - 0.5) - 90 for k in 1:ny]

"Cell-centre longitudes [°E] of an `nx`-column grid."
lons(nx::Integer) = [360 / nx * (i - 0.5) for i in 1:nx]

"""
    area_weights(nx, ny)

`cos(latitude)` cell weights, as the model's own global mean uses. A plain
`mean` over-weights the polar rows: on GREB's control Ts it reads 8.5-10.6 K
colder than the weighted mean.
"""
area_weights(nx::Integer, ny::Integer) = [cosd(φ) for _ in 1:nx, φ in lats(ny)]

"""
    by_position(n)

The `(year, month)` of each of `n` records that carry no times: January of year
1 onward. A `greb_model!` result brings its own times (`ctrl_time`,
`scnr_time`); this is for a bare record vector.
"""
by_position(n::Integer) = [(year=cld(i, 12), month=mod1(i, 12)) for i in 1:n]

"Record indices of each whole year in `times`, in order; a year without all 12 months is dropped."
function yearblocks(times::AbstractVector)
    blocks = [findall(t -> t.year == y, times) for y in unique(t.year for t in times)]
    filter(b -> length(b) == 12, blocks)
end

"""
    series(records, var; mask=nothing)

Area-weighted mean of `var`, one value per record: the global mean, or the mean
over the cells where `mask` (a lon x lat `Bool` matrix, see
[`band_mask`](@ref) and [`land_mask`](@ref)) is true.
"""
function series(records::AbstractVector, var::Symbol; mask::Union{Nothing,AbstractMatrix{Bool}}=nothing)
    w = area_weights(size(getfield(first(records), var))...)
    mask === nothing || (w = w .* mask)
    sum(w) > 0 || error("the mask selects no cells")
    [sum(getfield(r, var) .* w) / sum(w) for r in records]
end

"""
    to_absolute(anomaly, control)

A scenario that is an anomaly, made absolute again: the control's last 12 months
are added back by calendar month. Needed to difference a scenario against its
control with [`plot_difference`](@ref).
"""
function to_absolute(anomaly::AbstractVector, control::AbstractVector)
    length(control) >= 12 || error("the control needs at least 12 months, it has $(length(control))")
    clim = control[end-11:end]
    return [map(.+, r, clim[mod1(i, 12)]) for (i, r) in enumerate(anomaly)]
end

"""
    zonal_mean(records, var; month=:mean)

The mean over longitude of one grid ([`field`](@ref)), one value per latitude.
"""
zonal_mean(records::AbstractVector, var::Symbol; month=:mean) = vec(mean(field(records, var; month=month); dims=1))

"""
    annual(s, times=by_position(length(s)))

Calendar-year means of a monthly series, one per whole year. `times` holds the
`(year, month)` of each value.
"""
annual(s::AbstractVector, times::AbstractVector=by_position(length(s))) = [mean(s[b]) for b in yearblocks(times)]

"""
    seasonal_cycle(records, var, times=by_position(length(records)))

Each calendar month's global mean, averaged over the whole years present.
"""
function seasonal_cycle(records::AbstractVector, var::Symbol, times::AbstractVector=by_position(length(records)))
    s = series(records, var)
    whole = reduce(vcat, yearblocks(times); init=Int[])
    isempty(whole) ? Float64[] : [mean(s[i] for i in whole if times[i].month == m) for m in 1:12]
end

"""
    field(records, var; month=:mean)

One grid: month index `month`, the last month (`:last`), or the time-mean (`:mean`).
"""
function field(records::AbstractVector, var::Symbol; month=:mean)
    month === :mean && return mean(getfield(r, var) for r in records)
    month === :last && return getfield(records[end], var)
    getfield(records[month], var)
end

"""
    hovmoller(records, var) -> (months, latitudes, matrix)

Zonal mean per month: `matrix[t, j]` is month `t` at latitude `j`.
"""
function hovmoller(records::AbstractVector, var::Symbol)
    m = permutedims(reduce(hcat, [vec(mean(getfield(r, var); dims=1)) for r in records]))
    (1:length(records), lats(size(m, 2)), m)
end

"""
    map_frames(records, var; step=:month, times=by_position(length(records))) -> Vector{Matrix{Float64}}

One grid per month, or per whole calendar year with `step=:year`.
"""
function map_frames(records::AbstractVector, var::Symbol; step::Symbol=:month,
                    times::AbstractVector=by_position(length(records)))
    frames = [Float64.(getfield(r, var)) for r in records]
    step === :month && return frames
    step === :year && return [mean(frames[b]) for b in yearblocks(times)]
    error("step must be :month or :year, got :$step")
end
