# Reading the plots

## Absolute values and anomalies

The control run is drawn in absolute values on the viridis scale. A scenario that is an anomaly
(`scnr_anomaly` is true in the result) is the change against the control and is drawn on a scale
centred at zero: red is warmer or more, blue colder or less, white no change. A scenario that is not an
anomaly, such as an orbital run, is drawn in absolute values. A result without `scnr_anomaly`, and a bare
record vector, count as absolute: the plots do not guess.

Panel titles say which: `scenario anomaly - Surface temperature`.

## Colour scales

| Plot | Scale |
|:-----|:------|
| `plot_map`, `plot_hovmoller` | Set by the data of each panel; the control and the scenario do not share a scale |
| `plot_difference` | Symmetric about zero, set by the largest difference |
| Animations | Fixed for all frames of a run: the same colour is the same value in every frame |
| Anomaly animations | Symmetric about zero at the 99th percentile of the absolute values, so that a few cells with large changes do not wash out the rest. Those cells saturate |

Compare panels by their colour bar, not by colour alone.

## Means and weights

A global mean weights each cell by `cos(latitude)`, because cells near the poles cover less area. A plain
average over the grid counts every row equally and reads about 8.5 to 10.6 K too cold for the control surface
temperature. All the plots and `series` use the weighted mean; masks keep the weights inside the region.

## What a 96x48 grid can show

| Feature | What to expect |
|:--------|:---------------|
| Cell size | 3.75 by 3.75 degrees. Features smaller than a few cells, such as a lake or a mountain range, are not resolved |
| Coastlines | The outline follows cell edges: land is where the model's topography is above 0 m, so small islands and straits are missing or wide |
| Poles | The first and last rows are bands near the poles. Zonal means there average over short circles |
| Smoothness | Maps are drawn cell by cell, not smoothed. Sharp edges in a map are cell edges |

## Where a signal comes from

| Question | Plot to use |
|:---------|:------------|
| Where did it change? | `plot_map` of the scenario anomaly, or `plot_difference` of two runs |
| How much, at which latitude? | `plot_zonal` |
| When in the year? | `plot_seasonal`, or `plot_hovmoller` for latitude and season together |
| How fast, and has it settled? | `plot_timeseries` with `annual = true`: a flat line at the end means the run has reached a new balance, a rising line that it has not |
| Which process matters? | Run the model with the process switched off and use `plot_compare` or `plot_members` (see `examples/feedback_deconstruction.jl`) |

## Limits to keep in mind

- The control run drifts little by design: the model is flux-corrected against the climatology during the
  spin-up. A flat control line is expected, not a sign of an error.
- A short scenario shows the fast response. After 15 years the warming is still rising (see the annual time
  series); `plot_timeseries` shows how much is left.
- Each member of an ensemble with different physics needs its own spin-up before its output means anything.
- The viz plots show model output; they do not test it. A plausible-looking map can still be wrong, so the
  run checks of the package (`RangeCheck`) come first.
