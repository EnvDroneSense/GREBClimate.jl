# Experiments

The scripts in `tools/experiments/` run the model to answer a question. They need the local dataset, read it with `greb_data_dir(; allow_download=false)` and never download it. Run each from the repository root.

| Script | Answers | Typical run time |
|:-------|:--------|:-----------------|
| `sensitivity.jl` | What is the effective climate sensitivity of an abrupt-CO2 preset? | Minutes for 100 years |
| `co2_sweep.jl` | How does the last-year climate change with the CO2 level? | Minutes for the default 5 levels of 30 years; hours for a full curve |
| `helpers.jl` | The module `ExperimentTools` that both scripts use | Not run |

## Climate sensitivity

```bash
julia --project=. -t 2,0 tools/experiments/sensitivity.jl [preset] [scenario years] [spin-up years] [control years]
```

The defaults are `co2_double`, 100 scenario years, a spin-up of 3 years and a control of 5 years, on the NCEP climatology. The script regresses the change in the top-of-atmosphere net flux (`sw - olr`) on the change in surface temperature over the annual global means (the Gregory method). The preset must give its scenario as an anomaly. A short run:

```text
co2_double, 20 scenario years, spin-up 2, control 3
  forcing                   4.95 W/m2
  feedback                 -1.81 W/m2/K
  effective sensitivity     2.73 K
  r2                       0.999
  warming reached           2.50 K
  net flux remaining        0.39 W/m2
```

| Line | Meaning |
|:-----|:--------|
| `forcing` | The net flux at zero warming, from the fit |
| `feedback` | The slope of the fit |
| `effective sensitivity` | The warming at which the fit reaches zero net flux: `-forcing / feedback` |
| `warming reached` | The last annual `Ts` change. It can differ from the sensitivity because `sw - olr` is not the whole budget of the anomaly |
| `net flux remaining` | The last annual net-flux change |

The 100-year default for `co2_double` gives 5.06 W/m2, -1.88 W/m2/K and 2.69 K, with 3.00 K reached.

## CO2 sweep

```bash
julia --project=. -t 4,0 tools/experiments/co2_sweep.jl [levels] [scenario years] [spin-up years] [control years]
```

The script holds CO2 at `levels` values (default 5), log-spaced between 170 ppm (`:co2_half`) and 3400 ppm (`:co2_10x`), after a control at 280 ppm. It needs at least 12 scenario years, since it takes a last-year mean. It prints the anomaly of the last scenario year and writes it to `tools/experiments/out/co2_sweep.csv` (git-ignored). A short run:

```text
 CO2 [ppm]    Ts [K]       ice    precip
     170.0     -1.27    0.0113    -0.176
     760.0      2.95   -0.0228     0.484
    3400.0      8.22   -0.0492      1.66
```

| Column | Unit |
|:-------|:-----|
| `Ts` | K, global mean |
| `ice` | Fraction of the grid, global mean |
| `precip` | mm/day, global mean |

The levels run side by side, which is why it wants several threads. Each running level holds 370 MB of fields, so check the free memory before you raise `-t`.

## ExperimentTools

`helpers.jl` defines the module `ExperimentTools`. Its functions are pure: none runs the model, and every mean is area-weighted. Include it with `include("tools/experiments/helpers.jl"); using .ExperimentTools: summarize, gregory`.

| Function | Returns |
|:---------|:--------|
| `annual_means(records, times, f)` | The years that have all 12 months, and the annual mean of a field name or a function of a record |
| `summarize(result)` | A few kilobytes: annual global means of `Ts`, `precip`, `ice` and the net flux for control and scenario, and the first record outside the allowed range |
| `in_range_report(records, times, anomaly)` | `nothing`, or the year, month, field and value of the first record outside the range |
| `gregory(dTs, dnet)` | `forcing`, `feedback`, `sensitivity`, `r2`, `reached` and `remaining` |
| `area_mean`, `zonal_mean` | The area-weighted mean of a map (optionally under a mask), and its zonal mean |
| `polar_amplification(map)` | The mean poleward of a latitude (default 60), north and south, divided by the global mean |
| `land_ocean(map, z_topo)` | The area-weighted means over land and over ocean |
| `pattern_correlation(a, b)` | The area-weighted correlation of two maps |
| `sine_fit(series, period)` | The amplitude, the lag in samples and the mean of a fitted sine |

`test/test_experiment_tools.jl` tests the helpers; run it by name with `Pkg.test(test_args=["experiment_tools"])` after changing them.
