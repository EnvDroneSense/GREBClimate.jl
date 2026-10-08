# The per-step observer of greb_model! and the budget check built on it.

function _budget_run(; steps, setup! = (Ts, Ta, q, To, fields) -> nothing, config = preset(:full_model))
    fields = synthetic_fields()
    r = resolve(config)
    ini = quiet(() -> init_model!(r, fields))
    fields.Ts_flux_correction .= 0.25f0
    fields.To_flux_correction .= 0.01f0
    fields.q_flux_correction .= 1.0f-7
    Ts, Ta, q, To = copy(ini.Ts_ini), copy(ini.Ta_ini), copy(ini.q_ini), copy(ini.To_ini)
    setup!(Ts, Ta, q, To, fields)
    check = GREBClimate.BudgetCheck()
    state, ws, acc, ts = ModelState(), ModelWorkspace(), MonthlyAccumulator(), TimeState(1, 1)
    quiet() do
        for step in 1:steps
            time_loop!(control_time(step), ini.CO2_ctrl, Ts, Ta, q, To, MonthlyRecord[], GREBClimate.RecordTime[],
                       fields, state, ws, acc, ts, r; observer = check)
        end
    end
    return check
end

@testset "greb_model! calls the observer before and after every step's update" begin
    calls = Tuple{GREBClimate.Phase,Symbol}[]
    its = Int[]
    times = ModelTime[]
    observer = (point, view) -> (push!(calls, (view.phase, point)); push!(its, view.step); push!(times, view.time))
    quiet() do
        greb_model!(RunSpec(ctrl = 1, scnr = 1), preset(:full_model; corrections = NoCorrections());
                    fields = synthetic_fields(), allow_uninitialized = true, observer = observer)
    end
    step = [:after_tendencies, :after_step]
    @test calls == [(phase, point) for phase in (GREBClimate.control, GREBClimate.scenario) for _ in 1:N for point in step]
    @test its == [it for _ in 1:2 for it in 1:N for _ in 1:2]
    # The view's time is the value the step ran at; the scenario starts in its own year
    @test first(times) === ModelTime(GREBClimate.control, 1970, 1)
    @test last(times) === ModelTime(GREBClimate.scenario, 1950, N)
end

@testset "the observer sees the state before the update, then after it" begin
    fields = synthetic_fields()
    r = resolve(preset(:full_model))
    ini = quiet(() -> init_model!(r, fields))
    Ts = copy(ini.Ts_ini)
    seen = Dict{Symbol,Matrix{Float32}}()
    quiet() do
        time_loop!(control_time(), ini.CO2_ctrl, Ts, copy(ini.Ta_ini), copy(ini.q_ini), copy(ini.To_ini),
                   MonthlyRecord[], GREBClimate.RecordTime[], fields, ModelState(), ModelWorkspace(), MonthlyAccumulator(),
                   TimeState(1, 1), r; observer = (point, view) -> (seen[point] = copy(view.Ts)))
    end
    @test seen[:after_tendencies] == ini.Ts_ini
    @test seen[:after_step] == Ts
    @test Ts != ini.Ts_ini
end

@testset "BudgetCheck: each store changes by the sum of its flows" begin
    # A few ulp: 1e-4 K at 280 K, 1e-8 kg/kg at 0.005. A flux error of 1 W/m2
    # over the ocean moves Ts by 2e-4 K in one step.
    check = _budget_run(steps = 20)
    @test check.steps == 20
    @test check.surface < 1e-4
    @test check.atmosphere < 1e-4
    @test check.ocean < 1e-4
    @test check.water < 1e-8
    @test check.floor_Ts == 0
    @test check.floor_Ta == 0
    @test check.humidity_high == 0
    @test check.rain_limit == 0
end

@testset "BudgetCheck counts the cells a limiter held" begin
    # Far enough below the floor that one step cannot lift them back over it
    cold = _budget_run(steps = 1, setup! = function (Ts, Ta, q, To, fields)
        Ts[60, 4] = GREBClimate.min_T_K - 0.5f0
        fields.cap_surf[60, 4] = 1.0f12
    end)
    @test (cold.floor_Ts, cold.floor_Ta) == (1, 0)
    @test cold.surface < 1e-4          # the floor is part of the expected value
    cold = _budget_run(steps = 1, setup! = (Ts, Ta, q, To, fields) -> (Ta[5, 6] = -5000.0f0))
    @test (cold.floor_Ts, cold.floor_Ta) == (0, 1)

    wet = _budget_run(steps = 1,
        setup! = (Ts, Ta, q, To, fields) -> (fields.q_flux_correction[7, 8, 1] = 0.5f0; fields.q_flux_correction[9, 10, 1] = -0.5f0))
    # The fixture's coldest row rains out more than it holds on its own
    base = _budget_run(steps = 1)
    @test wet.humidity_high == base.humidity_high + 1
    @test wet.humidity_low == base.humidity_low + 1
    @test wet.water < 1e-8

    limited = _budget_run(steps = 1, config = preset(:full_model; hydrology = Hydrology(rain = :rh)),
        setup! = (Ts, Ta, q, To, fields) -> (fields.rain_limit .= -1.0f0))
    @test limited.rain_limit == X * Y
end

@testset "RangeCheck: the first step outside the physical range is recorded" begin
    G = GREBClimate
    field(v) = fill(Float32(v), X, Y)
    view(step; Ts = 280, Ta = 270, To = 285, q = 0.005) =
        (; phase = G.control, step, year = 1950, Ts = field(Ts), Ta = field(Ta), To = field(To), q = field(q))

    check = G.RangeCheck()
    check(:after_tendencies, view(1; Ts = 1e24))     # only the state after the step is checked
    @test check.steps == 0 && G.in_range(check)
    check(:after_step, view(1))
    check(:after_step, view(2; Ta = 265))
    @test check.steps == 2 && G.in_range(check)
    @test check.seen.Ta == (265.0f0, 270.0f0) && check.seen.Ts == (280.0f0, 280.0f0)

    # A cell that ran away, then a later one: the first is kept
    runaway = view(3)
    runaway.Ts[5, 7] = 1.0f24
    check(:after_step, runaway)
    check(:after_step, view(4; q = -1))
    @test !G.in_range(check)
    @test check.first == (phase = G.control, step = 3, year = 1950, field = :Ts, value = 1.0f24)
    @test check.seen.q[1] == -1.0f0

    # A NaN is outside; limits can be set
    nan = G.RangeCheck()
    nan(:after_step, view(1; To = NaN))
    @test !G.in_range(nan) && nan.first.field === :To
    narrow = G.RangeCheck(Ts = (285, 300))
    narrow(:after_step, view(1))
    @test narrow.first.field === :Ts && narrow.first.value == 280.0f0

    # In a run it sees every step of the control and the scenario
    seen = G.RangeCheck(Ts = (-Inf, Inf), Ta = (-Inf, Inf), To = (-Inf, Inf), q = (-Inf, Inf))
    quiet() do
        greb_model!(RunSpec(ctrl = 1, scnr = 1), preset(:full_model; corrections = NoCorrections());
                    fields = synthetic_fields(), allow_uninitialized = true, observer = seen)
    end
    @test seen.steps == 2N
end
