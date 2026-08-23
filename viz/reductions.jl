# Pure reductions: model output in, numbers out. No plotting, no model import -
# everything here works on plain NamedTuples of matrices, and the grid size is
# read from the data rather than assumed. That keeps the tool independent of
# how (or whether) GREBClimate was loaded.

# ── Grid geometry ──────────────────────────────────────────────────────────

"Cell-centre latitudes [deg] for a grid `ny` rows tall (matches `lat_grid`)."
lats(ny::Integer) = [180 / ny * k - (180 / ny) / 2 - 90 for k in 1:ny]

"Cell-centre longitudes [deg east] for a grid `nx` columns wide."
lons(nx::Integer) = [360 / nx * i - (360 / nx) / 2 for i in 1:nx]

"""
    area_weights(nx, ny) -> Matrix{Float64}

`cos(latitude)` cell weights. A lat-lon grid has as many cells at 88° as at the
equator, so a plain `mean` over-weights the poles badly - measured at ~10 K on
GREB's global-mean Ts. This is the weighting Fortran `gmean()` applies
(greb.model.mscm.f90:1497-1513), and what the model's own `diagnostics!` uses.
"""
area_weights(nx::Integer, ny::Integer) = [cosd(φ) for _ in 1:nx, φ in lats(ny)]

"Grid size of `records`, read from the first record's `var` field."
function griddims(records, var::Symbol)
    isempty(records) && error("no records")
    size(getfield(first(records), var))
end

# ── Default regions ────────────────────────────────────────────────────────
# Each is one `register_region!` call. Add your own the same way; the notebook
# picks it up without being edited.

_latgrid(nx, ny) = [φ for _ in 1:nx, φ in lats(ny)]
_needfields(f, name) = f === nothing ?
    error("region :$name needs `fields=` (it reads z_topo)") : f

register_region!(:global,    "Global",              (nx, ny, _) -> trues(nx, ny))
register_region!(:nh,        "N. hemisphere",       (nx, ny, _) -> _latgrid(nx, ny) .> 0)
register_region!(:sh,        "S. hemisphere",       (nx, ny, _) -> _latgrid(nx, ny) .< 0)
register_region!(:tropics,   "Tropics (±23.5°)",    (nx, ny, _) -> abs.(_latgrid(nx, ny)) .<= 23.5)
register_region!(:extratrop, "Extratropics",        (nx, ny, _) -> abs.(_latgrid(nx, ny)) .> 23.5)
register_region!(:arctic,    "Arctic (>66.5°N)",    (nx, ny, _) -> _latgrid(nx, ny) .>= 66.5)
register_region!(:antarctic, "Antarctic (>66.5°S)", (nx, ny, _) -> _latgrid(nx, ny) .<= -66.5)
register_region!(:land,  "Land only",  (nx, ny, f) -> _needfields(f, :land).z_topo .> 0)
register_region!(:ocean, "Ocean only", (nx, ny, f) -> _needfields(f, :ocean).z_topo .< 0)

# ── Time series ────────────────────────────────────────────────────────────

"""
    series(records, var; region=:global, weighted=true, fields=nothing)

Collapse each monthly record to one number, giving a value per month.

Area-weighted by default. `weighted=false` reproduces a plain `mean`, which is
not a physically meaningful average - it exists only for comparing against code
that does it the old way.
"""
function series(records::AbstractVector, var::Symbol;
                region::Symbol=:global, weighted::Bool=true, fields=nothing)
    isempty(records) && return Float64[]
    nx, ny = griddims(records, var)
    w = (weighted ? area_weights(nx, ny) : ones(nx, ny)) .*
        region_mask(region, nx, ny; fields=fields)
    tot = sum(w)
    tot == 0 && error("region :$region selected no cells")
    [sum(getfield(r, var) .* w) / tot for r in records]
end

"""
    annual(s)

Average a monthly series into calendar years. A trailing partial year is
dropped rather than reported as a whole one.
"""
annual(s::AbstractVector) = [sum(@view s[i:i+11]) / 12 for i in 1:12:(length(s) - 11)]

# ── Spatial ────────────────────────────────────────────────────────────────

"""
    field(records, var; month=:mean)

One grid. `month` is an index, `:last`, or `:mean` (time-mean over every
record, usually what you want for a run longer than a year).
"""
function field(records::AbstractVector, var::Symbol; month=:mean)
    isempty(records) && error("no records")
    month === :last && return getfield(records[end], var)
    month === :mean && return sum(getfield(r, var) for r in records) ./ length(records)
    getfield(records[clamp(Int(month), 1, length(records))], var)
end

"""
    zonal_mean(records, var; month=:mean)

Longitude averaged away, leaving one value per latitude. No weighting needed:
cells in the same row have identical area.
"""
function zonal_mean(records::AbstractVector, var::Symbol; month=:mean)
    g = field(records, var; month=month)
    nx, ny = size(g)
    [sum(@view g[:, j]) / nx for j in 1:ny]
end

"""
    seasonal_cycle(records, var; region=:global, fields=nothing)

The 12-month climatology: each calendar month averaged over every year present.
Whole years only, so an incomplete final year cannot skew one month.
"""
function seasonal_cycle(records::AbstractVector, var::Symbol;
                        region::Symbol=:global, fields=nothing)
    s = series(records, var; region=region, fields=fields)
    nyr = length(s) ÷ 12
    nyr == 0 && return Float64[]
    [sum(s[m + 12k] for k in 0:(nyr - 1)) / nyr for m in 1:12]
end

"""
    hovmoller(records, var) -> (times, latitudes, matrix)

Zonal mean against both time and latitude, so a signal can be watched migrating
across latitudes - polar amplification, seasonal marches.
"""
function hovmoller(records::AbstractVector, var::Symbol)
    isempty(records) && error("no records")
    nx, ny = griddims(records, var)
    n = length(records)
    m = Matrix{Float64}(undef, n, ny)
    for (t, r) in enumerate(records)
        g = getfield(r, var)
        for j in 1:ny
            m[t, j] = sum(@view g[:, j]) / nx
        end
    end
    (1:n, lats(ny), m)
end

# ── Run helpers ────────────────────────────────────────────────────────────

"""
    difference(a, b, var; month=:mean)

`a - b` as a grid, for difference maps.
"""
difference(a::AbstractVector, b::AbstractVector, var::Symbol; month=:mean) =
    field(a, var; month=month) .- field(b, var; month=month)

# `records_of` / `label_of` / `kind_of` live in runs.jl.
