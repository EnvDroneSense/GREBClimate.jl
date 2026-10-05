# The calendar: 365 days, 2 steps a day, 730 steps a year, no leap years.
using Test
using GREBClimate
using GREBClimate: step_of_year, day_of_year, month_of_day, month_of_step,
    is_day_end, is_month_end, is_year_end, steps_in_month, first_step_of,
    decimal_year, MONTHS_PER_YEAR

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
        @test MONTHS_PER_YEAR == 12
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
