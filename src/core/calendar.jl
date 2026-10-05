# The model calendar: a 365-day year without leap years, stepped in 12-hour
# steps. Everything that turns a step number into a day, a month or a position
# in the year is here.

const ndays_yr = 365                            # days per year (no leap years)
const ndt_days = Int(round(24 * 3600 / Δt))     # time steps per day
"""
    nstep_yr

Time steps per year (`ndays_yr * ndt_days` = 730). Together with [`xdim`](@ref)
and [`ydim`](@ref) this fixes the shape of every field the model steps.

```jldoctest
julia> (xdim, ydim, nstep_yr)
(96, 48, 730)
```
"""
const nstep_yr = Int(ndays_yr * ndt_days)

const cjday_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]  # days per month
const jday_mon_cumsum = cumsum(cjday_mon)                           # last day of each month
const months_per_year = length(cjday_mon)

"""
    step_of_year(it) -> Int

Position of step `it` in the year, 1 to [`nstep_yr`](@ref). Steps are counted
from 1 at the first half of 1 January. It is also the index of the climatology
slice the model reads at that step.

```jldoctest
julia> step_of_year(1), step_of_year(730), step_of_year(731)
(1, 730, 1)
```
"""
step_of_year(it::Integer) = mod(it - 1, nstep_yr) + 1

"""
    day_of_year(it) -> Int

Day of the year of step `it`, 1 to 365.

```jldoctest
julia> day_of_year(1), day_of_year(2), day_of_year(3)
(1, 1, 2)
```
"""
day_of_year(it::Integer) = mod((it - 1) ÷ ndt_days, ndays_yr) + 1

"Month (1-12) that day `day` of the year (1-365) falls in."
function month_of_day(day::Integer)
    month = 1
    while day > jday_mon_cumsum[month]
        month += 1
    end
    return month
end

"""
    month_of_step(it) -> Int

Month (1-12) of step `it`.

```jldoctest
julia> month_of_step(62), month_of_step(63)
(1, 2)
```
"""
month_of_step(it::Integer) = month_of_day(day_of_year(it))

"True at the last step of a day."
is_day_end(it::Integer) = it % ndt_days == 0

"True at the last step of a month."
is_month_end(it::Integer) =
    is_day_end(it) && day_of_year(it) == jday_mon_cumsum[month_of_step(it)]

"True at the last step of a year."
is_year_end(it::Integer) = step_of_year(it) == nstep_yr

"Number of steps in `month` (1-12)."
steps_in_month(month::Integer) = cjday_mon[month] * ndt_days

"Step of the year of the first step of day `day` of `month`: `first_step_of(10, 1)` is 547."
function first_step_of(month::Integer, day::Integer)
    1 <= month <= months_per_year ||
        throw(ArgumentError("month must be 1 to $months_per_year, got $month"))
    1 <= day <= cjday_mon[month] ||
        throw(ArgumentError("month $month has days 1 to $(cjday_mon[month]), got $day"))
    days_before = month == 1 ? 0 : jday_mon_cumsum[month-1]
    return (days_before + day - 1) * ndt_days + 1
end

"""
    decimal_year(year, it) -> Float64

The time of step `it` in calendar `year` as a decimal year on the model
calendar: `year` at the first step of the year, rising by `1 / nstep_yr` per
step.

```jldoctest
julia> decimal_year(1991, 1), decimal_year(1991, 366)
(1991.0, 1991.5)
```
"""
decimal_year(year::Integer, it::Integer) = year + (step_of_year(it) - 1) / nstep_yr
