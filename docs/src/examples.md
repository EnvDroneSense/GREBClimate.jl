# Examples

Two scripts in `examples/` run the steps of the [Tutorial](tutorial.md) end to end. Run each from the repository root. Both read the dataset the way [`greb_data_dir`](@ref) does, and take an optional dataset directory as the last argument.

| Script | Question | Command |
|:-------|:---------|:--------|
| `examples/quickstart.jl` | What does one run look like, and what does a switched-off process change? | `julia --project=. -t 2,0 examples/quickstart.jl` |
| `examples/feedback_deconstruction.jl` | Which processes carry the warming of a doubled CO2? | `julia --project=. -t 4,0 examples/feedback_deconstruction.jl` |

## Quick start

The script runs `preset(:co2_double)` for 2 control and 10 scenario years, then runs it again with the ice-albedo feedback switched off. It prints the global-mean control temperature and the warming of each scenario year, which is an anomaly against the last control year.

```text
control: global-mean Ts 287.0 K

year   warming [K]   without ice albedo [K]
   1          0.62              0.51
   2          1.03              0.88
   3          1.21              1.06
   ...
  10          1.95              1.77
```

The listing is shortened. The code has five steps, one per section of the tutorial:

| Step | Code |
|:-----|:-----|
| 1. Input data | `dir = greb_data_dir(...)`, then `fields = load_climatology(dir; dataset=:ncep)` |
| 2. Experiment | `cfg = preset(:co2_double)` |
| 3. Run | `greb_model!(RunSpec(ctrl=2, scnr=10), cfg; jld2_dir=dir, fields=fields)` |
| 4. Result | `result.ctrl` and `result.scnr` hold one record per month |
| 5. A change | `preset(:co2_double; processes=(ice_albedo=false,))` runs the same experiment with one process off |

Pass `fields` to every run in a session so that the climatology is loaded once. With `--data-free` the script runs on empty fields and no dataset, which checks that the calls work; its numbers mean nothing.

## Feedback deconstruction

The script runs CO2 doubled once with the full model and once with each process switched off, side by side with [`run_ensemble`](@ref), and prints the warming of the last scenario year of each member. It writes the table to `examples/out/feedback_deconstruction.csv` (git-ignored).

```text
member                     warming [K]  vs full [K]
full model                        2.76          0.0
no ice albedo                     2.52        -0.24
uniform clouds                    2.58        -0.18
no moisture convergence           2.53        -0.23
no deep ocean                     2.99         0.23
no transport                      2.85         0.09
no water cycle                    1.58        -1.18
```

It takes about a minute with 4 threads. Members with different physics get their own spin-up. The differences do not add up to the full warming, because the processes interact. Each member holds about 370 MB of fields, so check the free memory before you raise the thread count.

## Further reading

| Page | Holds |
|:-----|:------|
| [Tutorial](tutorial.md) | The steps of a run, configuration switches, experiment keywords, results, ensembles |
| [Configuration](configuration.md) | Every `Config` part and its options |
| [Plots and notebook](viz/index.md) | The `viz/` toolbox that draws a result |
| [Experiments](dev/experiments.md) | The maintainer scripts that measure climate sensitivity and sweep the CO2 level |
