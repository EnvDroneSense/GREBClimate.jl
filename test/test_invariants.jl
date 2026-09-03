# The two properties the physics kernels rest on: allocations are a bounded
# constant regardless of grid size, and the return type is concrete. Both
# tables below are hand-written, so the first testset asserts they cover every
# kernel defined in the physics/circulation/tendencies files - add an eleventh
# kernel and that assertion fails rather than the kernel going unguarded.
#
# Scope is the physics kernels only. The other per-timestep functions -
# `output!`, `time_loop!`, `accumulate!`, `diagnostics!`, `forcing` - are not
# covered here. Package-level hygiene lives in test_aqua.jl.

# One fixture for both testsets. `fields` is an unloaded (all-zero) climatology
# on purpose - these tests measure allocations and inferred types, not physics,
# and an unloaded ClimateFields costs nothing to build.
function _kernel_fixture()
    X, Y = GREBClimate.xdim, GREBClimate.ydim
    fields = ClimateFields()
    cfg = quiet() do
        c = create_experiment_config(:full_model)
        init_model!(c, fields)
        c
    end
    return (
        fields = fields,
        cfg = cfg,
        state = ModelState(),
        ws = CirculationWorkspace(),
        ts = TimeState(1, 1),
        Ts = fill(290.0f0, X, Y),
        Ta = fill(280.0f0, X, Y),
        To = fill(285.0f0, X, Y),
        q = fill(0.006f0, X, Y),
        T1 = fill(280.0f0, X, Y),
        dX = zeros(Float32, X, Y),
    )
end

# Kernels write into pre-allocated `CirculationWorkspace` buffers and return a
# NamedTuple referencing them, so the only allocation is the tuple box - bytes,
# not scaling with xdim*ydim. A regression that allocates per grid cell
# overshoots by ~4 orders of magnitude, so scale matters, not exact counts.
@testset "physics kernels allocate a bounded constant" begin
    f = _kernel_fixture()

    # (name, thunk, byte budget). Budgets are the measured cost rounded up to
    # the next power of two. tendencies! is larger only because it returns a
    # 16-field tuple: with ws_a === ws_q === ws its two circulation! calls run
    # serially, so nothing here pays for @spawn.
    kernels = [
        ("SWradiation!", () -> SWradiation!(f.Ts, f.fields, f.state, f.ts, f.cfg, f.ws), 64),
        ("LWradiation!", () -> LWradiation!(f.Ts, f.Ta, f.q, 340.0f0, f.fields, f.ts, f.cfg, f.ws), 64),
        ("hydro!", () -> hydro!(f.Ts, f.q, f.fields, f.ts, f.cfg, f.ws), 64),
        ("seaice!", () -> seaice!(f.Ts, f.fields, f.ts, f.cfg), 0),
        ("deep_ocean!", () -> deep_ocean!(f.Ts, f.To, f.fields, f.ts, f.cfg, f.ws), 64),
        ("convergence!", () -> convergence!(f.T1, f.fields, f.ts, f.ws), 0),
        ("diffusion!", () -> diffusion!(f.T1, GREBClimate.z_air, f.fields, f.ws, f.ts), 0),
        ("advection!", () -> advection!(f.T1, GREBClimate.z_air, f.fields, f.ws, f.ts, f.cfg), 0),
        ("circulation!", () -> circulation!(f.T1, GREBClimate.z_air, f.dX, f.fields, f.ws, f.ts, f.cfg), 0),
        ("tendencies!", () -> tendencies!(340.0f0, f.Ts, f.Ta, f.To, f.q, f.fields, f.state, f.ws, f.ts, f.cfg), 256),
    ]

    # The table is hand-written; this is what stops a new kernel slipping past
    # it. "Kernel" means a `!`-named function defined in one of these files.
    kernel_files = ("physics/radiation.jl", "physics/hydrology.jl",
        "physics/ocean.jl", "circulation.jl", "tendencies.jl")
    defined_in_kernel_file(fn) = any(methods(fn)) do m
        path = replace(String(m.file), '\\' => '/')
        any(kf -> endswith(path, kf), kernel_files)
    end
    all_kernels = Set(String(n) for n in names(GREBClimate)
                      if endswith(String(n), "!") &&
                         getproperty(GREBClimate, n) isa Function &&
                         defined_in_kernel_file(getproperty(GREBClimate, n)))
    @test setdiff(all_kernels, Set(first.(kernels))) == Set{String}()

    for (name, thunk, budget) in kernels
        thunk()  # warm up: the first call pays for compilation
        # Best of two - GC bookkeeping can land on either call.
        allocated = min(@allocated(thunk()), @allocated(thunk()))
        @testset "$name <= $budget bytes (got $allocated)" begin
            @test allocated <= budget
        end
    end
end

# Type instability here is silent - the model just runs slower. A mutating
# kernel whose last expression is an `@.` broadcast returns Union{Nothing,Matrix}
# instead of Nothing; that is why seaice! ends in an explicit `return nothing`.
@testset "physics kernels have concrete return types" begin
    f = _kernel_fixture()

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

    for (kernel, argtypes) in signatures
        rts = Base.return_types(kernel, argtypes)
        @testset "$(nameof(kernel))" begin
            @test !isempty(rts)                      # the signature still matches a method
            @test all(isconcretetype, rts)
        end
    end

    # tendencies! matters most and is easiest to check on the real call:
    # @inferred throws unless the inferred type is concrete.
    @test (@inferred tendencies!(340.0f0, f.Ts, f.Ta, f.To, f.q, f.fields, f.state, f.ws, f.ts, f.cfg)) isa NamedTuple
end
