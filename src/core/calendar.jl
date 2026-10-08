# The model calendar: a 365-day year without leap years, stepped in 12-hour
# steps. Everything that turns a step number into a day, a month or a position
# in the year is here.

const days_per_year = 365                       # no leap years
const steps_per_day = Int(round(24 * 3600 / Δt))
"""
    nstep_yr

Time steps per year (`days_per_year * steps_per_day` = 730). Together with
[`xdim`](@ref) and [`ydim`](@ref) this fixes the shape of every field the model
steps.

```jldoctest
julia> (xdim, ydim, nstep_yr)
(96, 48, 730)
```
"""
const nstep_yr = Int(days_per_year * steps_per_day)

const days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
const last_day_of_month = cumsum(days_in_month)   # as a day of the year
const months_per_year = length(days_in_month)

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
day_of_year(it::Integer) = mod((it - 1) ÷ steps_per_day, days_per_year) + 1

"Month (1-12) that day `day` of the year (1-365) falls in."
function month_of_day(day::Integer)
    month = 1
    while day > last_day_of_month[month]
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
is_day_end(it::Integer) = it % steps_per_day == 0

"True at the last step of a month."
is_month_end(it::Integer) =
    is_day_end(it) && day_of_year(it) == last_day_of_month[month_of_step(it)]

"True at the last step of a year."
is_year_end(it::Integer) = step_of_year(it) == nstep_yr

"Number of steps in `month` (1-12)."
steps_in_month(month::Integer) = days_in_month[month] * steps_per_day

"Step of the year of the first step of day `day` of `month`: `first_step_of(10, 1)` is 547."
function first_step_of(month::Integer, day::Integer)
    1 <= month <= months_per_year ||
        throw(ArgumentError("month must be 1 to $months_per_year, got $month"))
    1 <= day <= days_in_month[month] ||
        throw(ArgumentError("month $month has days 1 to $(days_in_month[month]), got $day"))
    days_before = month == 1 ? 0 : last_day_of_month[month-1]
    return (days_before + day - 1) * steps_per_day + 1
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

"""
    Phase

The part of a run a step belongs to: `spinup`, `control` or `scenario`.
"""
@enum Phase::UInt8 spinup control scenario

"""
    ModelTime(phase, start_year, step)

Where a run is: the [`Phase`](@ref), the calendar year the phase started in,
and the step counted from 1 at the start of the phase. Everything else is
derived from these three.

```jldoctest
julia> t = ModelTime(GREBClimate.scenario, 1950, 731);

julia> GREBClimate.year(t), GREBClimate.month(t), step_of_year(t)
(1951, 1, 1)
```
"""
struct ModelTime
    phase::Phase
    start_year::Int
    step::Int
end

step_of_year(t::ModelTime) = step_of_year(t.step)
day_of_year(t::ModelTime) = day_of_year(t.step)
is_month_end(t::ModelTime) = is_month_end(t.step)
is_year_end(t::ModelTime) = is_year_end(t.step)

"Month (1-12) of `t`."
month(t::ModelTime) = month_of_step(t.step)

"Calendar year of `t`. Any whole number: 0 and negative years count on as the others do."
year(t::ModelTime) = t.start_year + (t.step - 1) ÷ nstep_yr

"""
    decimal_year(t::ModelTime) -> Float64

The time of `t` as a decimal year on the model calendar.
"""
decimal_year(t::ModelTime) = decimal_year(year(t), t.step)

"Index of the climatology slice the model reads at `t`."
data_slice(t::ModelTime) = step_of_year(t)

"Index of the slice before `slice`; the year repeats, so the first follows the last."
previous_slice(slice::Integer) = slice > 1 ? slice - 1 : nstep_yr

"Slice `slice` of a climatology `A` of size `(xdim, ydim, nstep_yr)`, as a view."
clim_slice(A::AbstractArray{<:Any,3}, slice::Integer) = @view A[:, :, slice]

struct PhaseSteps
    phase::Phase
    start_year::Int
    nsteps::Int
end

"""
    steps(phase, start_year, nyears)

Every [`ModelTime`](@ref) of a phase of `nyears` years that starts on 1 January
of `start_year`, in order.
"""
steps(phase::Phase, start_year::Integer, nyears::Integer) = PhaseSteps(phase, start_year, nyears * nstep_yr)

Base.length(s::PhaseSteps) = s.nsteps
Base.eltype(::Type{PhaseSteps}) = ModelTime
Base.iterate(s::PhaseSteps, step::Int=1) =
    step > s.nsteps ? nothing : (ModelTime(s.phase, s.start_year, step), step + 1)

"""
    RecordTime

The calendar `year` and `month` (1-12) of one [`MonthlyRecord`](@ref): a
`NamedTuple`.
"""
const RecordTime = @NamedTuple{year::Int, month::Int}
