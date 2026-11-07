# Building blocks

The plots are built from a few functions that you can use for plots of your own.
`reductions.jl` and `regions.jl` need no plotting package; the plots add Plots.jl on top.
A record is a `NamedTuple` of lon x lat matrices, one per month; the grid size is read
from the data.

## Numbers from records

| Function | Returns |
|:---------|:--------|
| `series(records, var; mask)` | The area-weighted mean of `var` per record; over a mask when given |
| `annual(s, times)` | The calendar-year means of a monthly series, whole years only |
| `seasonal_cycle(records, var, times)` | The mean of each calendar month over the whole years |
| `field(records, var; month)` | One grid: the time mean (`:mean`), the last month (`:last`) or month `n` |
| `zonal_mean(records, var; month)` | The mean over longitude, one value per latitude |
| `hovmoller(records, var)` | The zonal mean per month, as `(months, latitudes, matrix)` |
| `map_frames(records, var; step)` | One grid per month or per year |
| `to_absolute(anomaly, control)` | A scenario anomaly with the control's calendar-month means added back |
| `lats(ny)`, `lons(nx)`, `area_weights(nx, ny)` | Cell-centre coordinates and `cos(latitude)` weights |

The weights matter: a plain `mean` over-weights the polar rows (see [Reading the plots](interpreting.md)).

## Masks and outlines

A mask is a `Bool` matrix on the grid. The same mask selects cells for a mean and draws an outline on a map.

| Function | Returns |
|:---------|:--------|
| `land_mask(z_topo)` | The land cells (topography above 0 m) |
| `band_mask(nx, ny; south, north)` | The cells in a latitude band; `band_mask(grid; ...)` takes the size of a grid |
| `mask_outline(mask)` | The edges of a mask as `(lon, lat)` degree vectors, separated by `NaN`, for any line plot |
| `outline!(p, mask; color, lw)` | Draws that outline on a map |
| `coastlines!(p, fields)` | `outline!` of the land mask; what the maps use when given `fields` |

```julia
land = land_mask(fields.z_topo)
tropics = band_mask(land; south = -23.5, north = 23.5)
series(result.scnr, :Ts; mask = land .& tropics)       # tropical land warming per month

p = plot_map(result.ctrl; var = :Ts)                   # one run: one panel
outline!(p, tropics; color = :white)                   # mark a region on a map
```

Masks combine with `.&`, `.|` and `.!`. A mask of your own (a box of cells, a basin) works the same way.

## A plot of your own

Combine a reduction with Plots.jl. The helpers `geomap` (a heatmap on degree axes), `colours` (a zero-centred scale for anomalies) and `fieldinfo` (label and unit of a variable) are in `GREBViz` but not exported: use them as `GREBViz.geomap`.

```julia
using Plots
s = series(result.scnr, :Ts; mask = land)
y = annual(s, result.scnr_time)
plot(y; xlabel = "year", ylabel = "K", label = "land")

p = GREBViz.geomap(field(result.scnr, :precip); fields = fields, title = "Precipitation change")
```

## Not included

| Idea | Why not |
|:-----|:--------|
| Real coastlines (Natural Earth) | Adds a data dependency and shows detail the 3.75 degree grid cannot resolve; the cell outline shows what the model sees as land |
| Another backend (Makie, for hover values) | Would rewrite the plots; Plots.jl with GR is enough for static figures and the notebook |
| A Gregory plot, a land/ocean split plot, a polar projection | Not needed yet; `series(...; mask)` already gives the numbers for the first two |
| One shared set of reductions with the package and `tools/` | `viz/` has its own environment and does not import the model, which keeps its tests free of the dataset |
