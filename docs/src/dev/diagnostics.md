# Diagnostics

`tools/diagnostics/budget.jl` checks the model's bookkeeping. For one preset it runs a control and a scenario phase with an observer attached and reports, per phase, whether each store changes by the sum of its flows, how often a limiter held a cell, and the global-mean flows. It needs the local dataset.

```bash
julia --project=. -t 2,0 tools/diagnostics/budget.jl [preset] [years]
```

The defaults are `full_model` and 1 year of control and 1 year of scenario, on the stored flux corrections. The observer is `GREBClimate.BudgetCheck` (not exported), which recomputes each step's expected `Ts`, `Ta`, `To` and `q` from the step's flows, limiters included, and compares them with what the model stored.

```text
preset :full_model, 1 yr control + 1 yr scenario, 2 thread(s)

== control (730 steps) ==
largest gap between a store's change and its flows:
  surface 3.05e-05 K   air 3.05e-05 K   ocean 0.00e+00 K   humidity 3.73e-09 kg/kg
limiter cell-steps (of 3363840):
  Ts floor 0   Ta floor 0   humidity removal 113336   humidity gain 0   rain limit 0
global means:
  shortwave absorbed by the surface [W/m2]         +233.596
  longwave out to space (olr) [W/m2]               +225.443
  net flux into the surface [W/m2]                 -27.2487
  ...
```

The scenario block has the same layout. The listing is shortened here.

| Line | Read it as |
|:-----|:-----------|
| `largest gap ...` | The largest difference over all steps and cells between the stored and the expected value. A gap near single-precision rounding (about 1e-5 K) means the store closes. A larger gap means a flow is missing from the update or counted twice |
| `Ts floor`, `Ta floor` | Cell-steps where a temperature floor held the value |
| `humidity removal`, `humidity gain` | Cell-steps where the humidity limiter removed or added water. Removal is common in the default run; a large gain is not expected |
| `rain limit` | Cell-steps where rain was set to the rain limit (only with `rain = :rh`) |
| `global means` | Area-weighted mean of each flow per step, in the unit in brackets |

## When to run it

Run it after a change to the flux sums, the limiters or the order of the tendency assembly. Compare the gaps and the limiter counts before and after. A change that must not alter results is checked with the [bit-identity tool](validation.md) instead.

`test/test_budgets.jl` holds the automatic version of this check on a synthetic grid. The script is for the real dataset and the real presets.
