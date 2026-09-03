# The two properties the whole kernel design rests on: allocations are a
# bounded constant regardless of grid size, and every kernel's return type is
# concrete. Package-level hygiene lives in test_aqua.jl.

# Kernels write into pre-allocated `CirculationWorkspace` buffers and return a
# NamedTuple referencing them, so the only allocation is the tuple box - bytes,
# not scaling with xdim*ydim. A regression that allocates per grid cell
# overshoots by ~4 orders of magnitude, so scale matters, not exact counts.
@testset "physics kernels allocate a bounded constant" begin
    X, Y = GREBClimate.xdim, GREBClimate.ydim
    fields = ClimateFields()
    state = ModelState()
    ws = CirculationWorkspace()
    ts = TimeState(1, 1)
    cfg = quiet() do
        c = create_experiment_config(:full_model)
        init_model!(c, fields)
        c
    end

    Ts = fill(290.0f0, X, Y)
    Ta = fill(280.0f0, X, Y)
    To = fill(285.0f0, X, Y)
    q = fill(0.006f0, X, Y)
    T1 = fill(280.0f0, X, Y)
    dX = zeros(Float32, X, Y)

    # (name, thunk, byte budget). Budgets are the measured cost rounded up to
    # the next power of two; tendencies! is larger because it returns a
    # 16-field tuple and may @spawn its two circulation! calls.
    kernels = [
        ("SWradiation!", () -> SWradiation!(Ts, fields, state, ts, cfg, ws), 64),
        ("LWradiation!", () -> LWradiation!(Ts, Ta, q, 340.0f0, fields, ts, cfg, ws), 64),
        ("hydro!", () -> hydro!(Ts, q, fields, ts, cfg, ws), 64),
        ("seaice!", () -> seaice!(Ts, fields, ts, cfg), 0),
        ("deep_ocean!", () -> deep_ocean!(Ts, To, fields, ts, cfg, ws), 64),
        ("convergence!", () -> convergence!(T1, fields, ts, ws), 0),
        ("diffusion!", () -> diffusion!(T1, GREBClimate.z_air, fields, ws, ts), 0),
        ("advection!", () -> advection!(T1, GREBClimate.z_air, fields, ws, ts, cfg), 0),
        ("circulation!", () -> circulation!(T1, GREBClimate.z_air, dX, fields, ws, ts, cfg), 0),
        ("tendencies!", () -> tendencies!(340.0f0, Ts, Ta, To, q, fields, state, ws, ts, cfg), 4096),
    ]

    for (name, f, budget) in kernels
        f()  # warm up: the first call pays for compilation
        # Best of two - GC bookkeeping can land on either call.
        allocated = min(@allocated(f()), @allocated(f()))
        @testset "$name <= $budget bytes (got $allocated)" begin
            @test allocated <= budget
        end
    end
end

# Type instability here is silent - the model just runs slower. A mutating
# kernel whose last expression is an `@.` broadcast returns Union{Nothing,Matrix}
# instead of Nothing; that is why seaice! ends in an explicit `return nothing`.
@testset "physics kernels have concrete return types" begin
    X, Y = GREBClimate.xdim, GREBClimate.ydim
    fields = ClimateFields()
    state = ModelState()
    ws = CirculationWorkspace()
    ts = TimeState(1, 1)
    cfg = quiet() do
        c = create_experiment_config(:full_model)
        init_model!(c, fields)
        c
    end

    F32 = Matrix{Float32}
    signatures = [
        (SWradiation!, (F32, ClimateFields, ModelState, TimeState, PhysicsConfig, CirculationWorkspace)),
        (LWradiation!, (F32, F32, F32, Float32, ClimateFields, TimeState, PhysicsConfig, CirculationWorkspace)),
        (hydro!, (F32, F32, ClimateFields, TimeState, PhysicsConfig, CirculationWorkspace)),
        (seaice!, (F32, ClimateFields, TimeState, PhysicsConfig)),
        (deep_ocean!, (F32, F32, ClimateFields, TimeState, PhysicsConfig, CirculationWorkspace)),
        (convergence!, (F32, ClimateFields, TimeState, CirculationWorkspace)),
        (diffusion!, (F32, Float32, ClimateFields, CirculationWorkspace, TimeState)),
        (advection!, (F32, Float32, ClimateFields, CirculationWorkspace, TimeState, PhysicsConfig)),
        (circulation!, (F32, Float32, F32, ClimateFields, CirculationWorkspace, TimeState, PhysicsConfig)),
    ]

    for (f, argtypes) in signatures
        rts = Base.return_types(f, argtypes)
        @testset "$(nameof(f))" begin
            @test length(rts) == 1
            @test isconcretetype(only(rts))
        end
    end

    # tendencies! matters most and is easiest to check on the real call:
    # @inferred throws unless the inferred type is concrete.
    Ts = fill(290.0f0, X, Y)
    Ta = fill(280.0f0, X, Y)
    To = fill(285.0f0, X, Y)
    q = fill(0.006f0, X, Y)
    @test (@inferred tendencies!(340.0f0, Ts, Ta, To, q, fields, state, ws, ts, cfg)) isa NamedTuple
end
