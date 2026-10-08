"""
    monthly_climatology(records::Vector{MonthlyRecord})::Vector{MonthlyRecord}

Returns a 12-month climatology taken from the *final* year of `records`.
"""
function monthly_climatology(records::Vector{MonthlyRecord})::Vector{MonthlyRecord}
    isempty(records) && return MonthlyRecord[]

    nmonths = months_per_year
    n = length(records)
    final_year = records[max(1, n - nmonths + 1):n]

    clim = MonthlyRecord[]
    for mon in 1:nmonths
        push!(clim, mon <= length(final_year) ? final_year[mon] : final_year[1])
    end
    return clim
end

"""
    scenario_anomalies(scnr_records, ctrl_clim)::Vector{MonthlyRecord}

Subtracts the matching calendar month of `ctrl_clim` (from
[`monthly_climatology`](@ref)) from each record in `scnr_records`,
turning absolute monthly output into anomalies relative to the control run.
"""
function scenario_anomalies(scnr_records::Vector{MonthlyRecord}, ctrl_clim::Vector{MonthlyRecord})::Vector{MonthlyRecord}
    isempty(scnr_records) && return scnr_records
    isempty(ctrl_clim) && return scnr_records

    anom = MonthlyRecord[]

    for (idx, rec) in enumerate(scnr_records)
        mon = mod(idx - 1, months_per_year) + 1
        ref = ctrl_clim[mon]
        push!(anom, map((x, x_ref) -> x .- x_ref, rec, ref))
    end
    return anom
end

"""
    ice_climatology(ctrl_output::Vector{MonthlyRecord})

Returns `ctrl_output`'s `ice` field from the *final* year, as an
`(xdim, ydim, 12)` array.
"""
function ice_climatology(ctrl_output::Vector{MonthlyRecord})
    ice_months = zeros(Float32, xdim, ydim, months_per_year)
    isempty(ctrl_output) && return ice_months

    nmonths = months_per_year
    n = length(ctrl_output)
    final_year = ctrl_output[max(1, n - nmonths + 1):n]

    for (mon, rec) in enumerate(final_year)
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
