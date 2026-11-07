# Masks and outlines on the lon x lat grid. A mask is a `Bool` matrix; the same
# mask selects cells for `series(...; mask)` and draws an outline on a map.

"The land cells, where the topography `z_topo` is above 0 m."
land_mask(z_topo::AbstractMatrix) = z_topo .> 0

"""
    band_mask(nx, ny; south=-90, north=90)
    band_mask(grid; south=-90, north=90)

The cells whose centre latitude lies in `[south, north]`, for an `nx` by `ny`
grid or the size of `grid`. `band_mask(96, 48; south=66.5)` is the Arctic.
"""
band_mask(nx::Integer, ny::Integer; south::Real=-90, north::Real=90) = [south <= φ <= north for _ in 1:nx, φ in lats(ny)]
band_mask(grid::AbstractMatrix; kw...) = band_mask(size(grid)...; kw...)

"""
    mask_outline(mask) -> (lon, lat)

The edges between the `true` and `false` cells of `mask`, as degree coordinates
for one line plot. Segments are separated by `NaN`. The longitude edge wraps
around the globe.
"""
function mask_outline(mask::AbstractMatrix{Bool})
    nx, ny = size(mask)
    dx, dy = 360 / nx, 180 / ny
    xs, ys = Float64[], Float64[]
    seg!(x1, y1, x2, y2) = (append!(xs, (x1, x2, NaN)); append!(ys, (y1, y2, NaN)))
    for j in 1:ny, i in 1:nx

        y0, y1 = -90 + (j - 1) * dy, -90 + j * dy
        mask[i, j] != mask[mod1(i + 1, nx), j] && seg!(i * dx, y0, i * dx, y1)       # east edge, wraps
        j < ny && mask[i, j] != mask[i, j+1] && seg!((i - 1) * dx, y1, i * dx, y1) # north edge
    end
    return xs, ys
end

"""
    outline!(p, mask; color=:black, lw=0.7)

Draw [`mask_outline`](@ref) on the map `p`. Drawn as line segments, not a contour:
a contour joins the heatmap's colour scale and is clipped when fixed `clims`
exclude its level.
"""
function outline!(p, mask::AbstractMatrix{Bool}; color=:black, lw=0.7)
    xs, ys = mask_outline(mask)
    Plots.plot!(p, xs, ys; linecolor=color, linewidth=lw, label="", xlims=(0, 360), ylims=(-90, 90))
end

"""
    coastlines!(p, fields; color=:black, lw=0.7)

Outline the land cells (`fields.z_topo > 0`) on the map `p`. No-op when `fields`
is `nothing`.
"""
coastlines!(p, fields; kw...) = fields === nothing ? p : outline!(p, land_mask(fields.z_topo); kw...)
