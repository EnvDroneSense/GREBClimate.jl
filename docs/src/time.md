# Time and calendar

```@meta
DocTestSetup = :(using GREBClimate)
```

The model has no dates and no clock. One value, a [`ModelTime`](@ref), says
where a run is, and the day, the month, the year and the climatology slice are
worked out from it.

## The calendar

| Unit | Length |
|:-----|:-------|
| Step | 12 hours |
| Day | 2 steps |
| Month | 28 to 31 days; February always has 28 |
| Year | 365 days, 730 steps ([`nstep_yr`](@ref)); no leap years |

Step 1 is the first half of 1 January. Every input climatology is one
repeating year of 730 slices, and the model reads slice `step_of_year` at each
step. Real dates are converted to this calendar before they reach the model;
the model itself never converts a date.

```jldoctest time
julia> t = ModelTime(GREBClimate.scenario, 1950, 1000);

julia> step_of_year(t), day_of_year(t), GREBClimate.month(t), GREBClimate.year(t)
(270, 135, 5, 1951)

julia> decimal_year(t)
1951.368493150685
```

[`step_of_year`](@ref), [`day_of_year`](@ref), [`month_of_step`](@ref) and
[`decimal_year`](@ref) also take a plain step number.

## The phases of a run

A run has up to three phases. Each starts on 1 January and counts its steps
from 1.

| Phase | Start year | Notes |
|:------|:-----------|:------|
| `GREBClimate.spinup` | 1970 | The flux-correction spin-up. Its log line counts `spin-up year 1, 2, ...` |
| `GREBClimate.control` | 1970 | The year is a label; nothing in the control depends on it |
| `GREBClimate.scenario` | `Scenario(start_year = ...)`, 1950 by default | The year drives the CO₂ path and the solar cycle |

The scenario does not continue from the control: both start from the same
initial state. A start year can be any whole number, including 0 and negative
years.

[`greb_model!`](@ref) returns the time of every record next to the records:

```julia
result.scnr_time[1]    # (year = 1950, month = 1)
```

An observer reads the time of the step from `view.time`.

## Where the stored year sits in the real year

The climatologies carry no dates. By label, slice 1 is the first half of
1 January and slice 730 the second half of 31 December.

The stored insolation table can be checked against an analytic calculation for
the modern orbit with the vernal equinox on day 80. The two agree best when one
is shifted by 1.5 days against the other (5.5 W/m² rms, against 6.7 W/m²
without a shift). The stored year is therefore aligned with the real year to
within about two days, not to the step. An input with real dates, such as an
eruption date, should not rely on a finer alignment than that.

## Not built

The time value was designed so that the following need no change to the
physics. None of them exists yet.

| Feature | What it would take |
|:--------|:-------------------|
| Start a run on another day | The phase would begin at a step other than 1, and the output would need a rule for the partial first month |
| Continue a stopped run | The time part is one `ModelTime`. The model state, the monthly sums and the surface heat capacity would have to be saved with it |
| CO₂ or sunlight that changes within the year | `decimal_year(t)` is available; today every forcing reads the whole-number year, except [`SeasonalCO2`](@ref), which reads the step of the year |
| Random draws that repeat in a continued run | A draw can be keyed on the phase and the step |
| Another step length | The steps per day and the map from time to climatology slice are the two places to change; the 730 slices of the dataset are the binding limit |
