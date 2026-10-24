"""
    monthly_climatology(records::Vector{MonthlyRecord})::Vector{MonthlyRecord}

Returns a 12-month climatology taken from the *final* year of `records`.
"""
function monthly_climatology(records::Vector{MonthlyRecord})::Vector{MonthlyRecord}
    isempty(records) && return MonthlyRecord[]
    final_year = _final_year(records)
    return [mon <= length(final_year) ? final_year[mon] : final_year[1] for mon in 1:months_per_year]
end

# The last twelve records, or all of them if there are fewer
_final_year(records) = records[max(1, end - months_per_year + 1):end]

"""
    scenario_anomalies(scnr_records, scnr_times, ctrl_clim)::Vector{MonthlyRecord}

Subtracts the matching calendar month of `ctrl_clim` (from
[`monthly_climatology`](@ref)) from each record in `scnr_records`,
turning absolute monthly output into anomalies relative to the control run.
`scnr_times` holds the [`RecordTime`](@ref) of each record.
"""
function scenario_anomalies(scnr_records::Vector{MonthlyRecord}, scnr_times::Vector{RecordTime},
    ctrl_clim::Vector{MonthlyRecord})::Vector{MonthlyRecord}
    (isempty(scnr_records) || isempty(ctrl_clim)) && return scnr_records
    return [map((x, x_ref) -> x .- x_ref, rec, ctrl_clim[time.month])
            for (rec, time) in zip(scnr_records, scnr_times)]
end

"""
    ice_climatology(ctrl_output::Vector{MonthlyRecord})

Returns `ctrl_output`'s `ice` field from the *final* year, as an
`(xdim, ydim, 12)` array.
"""
function ice_climatology(ctrl_output::Vector{MonthlyRecord})
    ice_months = zeros(Float32, xdim, ydim, months_per_year)
    isempty(ctrl_output) && return ice_months

    for (mon, rec) in enumerate(_final_year(ctrl_output))
        @. ice_months[:, :, mon] = rec.ice   # note: field name is `ice`, not `ice_cover`
    end

    return ice_months
end

"""
    global_mean(field::AbstractMatrix)

Area-weighted global mean of a field on the model grid (`xdim` by `ydim`):
each latitude row counts by the cosine of its latitude. A plain mean over the
cells counts the polar rows too heavily; for the surface temperature it reads
8 to 11 K colder.

```jldoctest
julia> global_mean(fill(288.0f0, xdim, ydim))
288.0
```
"""
function global_mean(field::AbstractMatrix)
    size(field) == (xdim, ydim) ||
        throw(DimensionMismatch("global_mean needs a ($xdim, $ydim) field, got $(size(field))"))
    total = 0.0
    for j in 1:ydim
        row = zero(eltype(field))
        @turbo for i in 1:xdim
            row += field[i, j]
        end
        total += Float64(row) * dxlat_grid[j]
    end
    return total / (xdim * sum(Float64, dxlat_grid))
end
