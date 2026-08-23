# State structs: constructors, field shapes, accumulator reset, mask reset.

@testset "co2_part regional CO2 mask resets between runs (no leak)" begin
    fields = ClimateFields()
    quiet() do
        init_model!(PhysicsConfig(experiment = :regional_co2_nh), fields)
    end
    @test any(!=(1.0f0), fields.co2_part)  # regional run actually changed the mask

    quiet() do
        init_model!(create_experiment_config(:full_model), fields)
    end
    @test all(==(1.0f0), fields.co2_part)  # init_model! resets it back to full CO2
end

@testset "workspace/accumulator/record types construct correctly" begin
    ws = CirculationWorkspace()
    @test size(ws.dTa_crcl) == (GREBClimate.xdim, GREBClimate.ydim)
    @test eltype(ws.dTa_crcl) === Float32

    acc = MonthlyAccumulator()
    fill!(acc.Tmm, 42.0f0)
    fill!(acc.qsensmm, 42.0f0)
    GREBClimate.reset!(acc)
    @test all(iszero, acc.Tmm)
    @test all(iszero, acc.qsensmm)

    ts = TimeState(1, 1)
    @test ts.jday == 1
    @test ts.ityr == 1

    @test MonthlyRecord <: NamedTuple
    @test :Ts in fieldnames(MonthlyRecord)
    @test :precip in fieldnames(MonthlyRecord)
end

@testset "state constructors: every field gets the right shape and eltype" begin
    X, Y, N = GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr

    # ClimateFields: 2D grid fields, the (ydim, nstep_yr) solar table,
    # the Bool flag, and everything else 3D.
    cf = ClimateFields()
    cf_2d = (:z_topo, :glacier, :z_ocean, :cap_surf, :wz_air, :wz_vapor,
             :rain_limit, :co2_part)
    @test length(fieldnames(ClimateFields)) == 39
    for f in fieldnames(ClimateFields)
        v = getfield(cf, f)
        if f === :loaded
            @test v === false
        elseif f === :sw_solar
            @test size(v) == (Y, N) && eltype(v) === Float32
        elseif f in cf_2d
            @test size(v) == (X, Y) && eltype(v) === Float32
        else
            @test size(v) == (X, Y, N) && eltype(v) === Float32
        end
    end
    # co2_part is the one field that is not zero-initialised.
    @test all(isone, cf.co2_part)
    for f in fieldnames(ClimateFields)
        f in (:loaded, :co2_part) && continue
        @test all(iszero, getfield(cf, f))
    end

    # CirculationWorkspace: four length-xdim vectors, the rest (xdim, ydim).
    cw = CirculationWorkspace()
    cw_vec = (:T1h, :dTxh, :term_north, :term_south)
    @test length(fieldnames(CirculationWorkspace)) == 42
    for f in fieldnames(CirculationWorkspace)
        v = getfield(cw, f)
        @test eltype(v) === Float32
        @test size(v) == (f in cw_vec ? (X,) : (X, Y))
        @test all(iszero, v)
    end

    # MonthlyAccumulator: 13 accumulators, all (xdim, ydim). There is no
    # `count` field - output! divides by cjday_mon[mon] * ndt_days.
    ma = MonthlyAccumulator()
    @test length(fieldnames(MonthlyAccumulator)) == 13
    for f in fieldnames(MonthlyAccumulator)
        v = getfield(ma, f)
        @test size(v) == (X, Y) && eltype(v) === Float32 && all(iszero, v)
    end

    # The keyword form is the point of the change: field-to-value
    # association by name, not by ordinal position.
    @test ClimateFields(loaded = true).loaded === true
    @test all(isone, ClimateFields(loaded = true).co2_part)
    @test size(MonthlyAccumulator(Tmm = zeros(Float32, 2, 2)).Tmm) == (2, 2)
    @test size(CirculationWorkspace(T1h = zeros(Float32, 3)).T1h) == (3,)
end
