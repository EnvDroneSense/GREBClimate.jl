# The calendar: 365 days, 2 steps a day, 730 steps a year, no leap years.
using Test
using GREBClimate
using GREBClimate: step_of_year, day_of_year, month_of_day, month_of_step,
    is_day_end, is_month_end, is_year_end, steps_in_month, first_step_of,
    decimal_year, months_per_year, spinup, control, scenario, steps, year, month, data_slice,
    previous_slice, clim_slice

@testset "calendar" begin
    @testset "step and day of the year" begin
        @test (step_of_year(1), day_of_year(1)) == (1, 1)
        @test (step_of_year(2), day_of_year(2)) == (2, 1)
        @test (step_of_year(730), day_of_year(730)) == (730, 365)
        @test (step_of_year(731), day_of_year(731)) == (1, 1)
        @test (step_of_year(1000), day_of_year(1000)) == (270, 135)
    end

    @testset "far past 200 years" begin
        @test (step_of_year(146_001), day_of_year(146_001)) == (1, 1)
        @test (step_of_year(73_000_000), day_of_year(73_000_000)) == (730, 365)
        @test (step_of_year(73_000_001), day_of_year(73_000_001)) == (1, 1)
    end

    @testset "months" begin
        @test months_per_year == 12
        @test month_of_day.([1, 31, 32, 59, 60, 90, 91, 365]) == [1, 1, 2, 2, 3, 3, 4, 12]
        @test month_of_step(62) == 1
        @test month_of_step(63) == 2
        @test steps_in_month.(1:12) == 2 .* [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        @test sum(steps_in_month, 1:12) == nstep_yr
    end

    @testset "ends of day, month and year" begin
        @test !is_day_end(1) && is_day_end(2)
        @test !is_month_end(61) && is_month_end(62) && !is_month_end(63)
        ends = filter(is_month_end, 1:(3 * nstep_yr))
        @test length(ends) == 36
        @test month_of_step.(ends) == repeat(1:12, 3)
        @test is_year_end(730) && !is_year_end(729) && is_year_end(1460)
    end

    @testset "dates to steps" begin
        @test first_step_of(1, 1) == 1
        @test first_step_of(4, 1) == 181
        @test first_step_of(10, 1) == 547
        @test first_step_of(12, 31) == 729
        @test_throws ArgumentError first_step_of(13, 1)
        @test_throws ArgumentError first_step_of(2, 29)
    end

    @testset "decimal year" begin
        @test decimal_year(1991, 1) == 1991.0
        @test decimal_year(1991, 730) == 1991 + 729 / 730
        @test decimal_year(1991, 731) == 1991.0   # the caller advances the year
        @test decimal_year(-5, 366) == -5 + 365 / 730
    end
end

@testset "model time" begin
    @testset "a plain value, everything derived" begin
        @test isbitstype(ModelTime)
        t = ModelTime(scenario, 1950, 1000)
        @test (step_of_year(t), day_of_year(t), month(t), year(t)) == (270, 135, 5, 1951)
        @test data_slice(t) == step_of_year(t)
        @test decimal_year(t) == 1951 + 269 / 730
        @test all(s -> month(ModelTime(control, 1970, s)) == month_of_step(s), 1:(2 * nstep_yr))
        @test all(s -> is_month_end(ModelTime(control, 1970, s)) == is_month_end(s), 1:(2 * nstep_yr))
        @test is_year_end(ModelTime(control, 1970, 730)) && !is_year_end(ModelTime(control, 1970, 731))
    end

    @testset "years at and below zero" begin
        for start in (1970, 1, 0, -231_000)
            counter = start
            for s in 1:(3 * nstep_yr)
                @test year(ModelTime(scenario, start, s)) == counter
                is_year_end(s) && (counter += 1)
            end
        end
        @test decimal_year(ModelTime(scenario, -5, 366)) == -5 + 365 / 730
    end

    @testset "the steps of a phase" begin
        phase = steps(scenario, 1950, 2)
        @test length(phase) == 2 * nstep_yr && eltype(phase) === ModelTime
        all_steps = collect(phase)
        @test first(all_steps) === ModelTime(scenario, 1950, 1)
        @test last(all_steps) === ModelTime(scenario, 1950, 2 * nstep_yr)
        @test year(last(all_steps)) == 1951
        # Walking a phase allocates nothing per step
        walk(p) = (n = 0; for t in p; n += year(t); end; n)
        walk(steps(scenario, 1950, 1))
        @test @allocated(walk(steps(scenario, 1950, 100))) == @allocated(walk(steps(scenario, 1950, 1)))
    end

    @testset "from time to data" begin
        @test data_slice(ModelTime(scenario, 1950, 731)) == 1
        @test previous_slice(2) == 1 && previous_slice(nstep_yr) == nstep_yr - 1
        @test previous_slice(1) == nstep_yr   # the year repeats
        A = reshape(collect(1:24), 2, 3, 4)
        @test clim_slice(A, 3) == A[:, :, 3]
        clim_slice(A, 3)[1, 1] = -1
        @test A[1, 1, 3] == -1                # a view, not a copy
    end

    @testset "empty phase" begin
        @test isempty(collect(steps(spinup, 1970, 0)))
        @test length(steps(control, 1970, 0)) == 0
    end
end
