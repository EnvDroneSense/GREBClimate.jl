# State structs: field shapes, accumulator reset, mask reset, derived fields.

@testset "init_model! clears the CO2 mask an earlier run left in the same fields" begin
    fields = ClimateFields()
    GREBClimate.apply_co2_mask!(LatitudeMask(:nh), fields)
    @test any(!=(1.0f0), fields.co2_part)  # what a regional scenario leaves behind

    # Full CO2 again, also for a regional preset: its mask is for the scenario
    for p in (:full_model, :regional_co2_nh)
        fields.co2_part[1, 1] = 0.5f0
        quiet() do
            init_model!(resolve(preset(p)), fields)
        end
        @test all(==(1.0f0), fields.co2_part)
    end
end

@testset "reset! zeroes every accumulator field" begin
    acc = MonthlyAccumulator()
    foreach(f -> fill!(getfield(acc, f), 42.0f0), fieldnames(MonthlyAccumulator))
    GREBClimate.reset!(acc)
    @test all(f -> all(iszero, getfield(acc, f)), fieldnames(MonthlyAccumulator))
end

@testset "state constructors: every field gets the right shape and eltype" begin
    X, Y, N = GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr

    # ClimateFields: 2D grid fields, the (ydim, nstep_yr) solar table,
    # the Bool flag, the two anomaly sources, and everything else 3D.
    cf = ClimateFields()
    cf_2d = (:z_topo, :glacier, :z_ocean, :cap_surf, :wz_air, :wz_vapor,
             :rain_limit, :co2_part)
    for f in fieldnames(ClimateFields)
        v = getfield(cf, f)
        if f === :loaded
            @test v === false
        elseif f === :anom_cc_source
            @test v == ""
        elseif f === :anom_enso_source
            @test v == ("", :none)
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
        f in (:loaded, :co2_part, :anom_cc_source, :anom_enso_source) && continue
        @test all(iszero, getfield(cf, f))
    end

    # CirculationWorkspace: four vectors, the rest matrices. The zonal-stencil
    # buffers carry longitude ghost cells, so their first dimension is `xghost`.
    cw = CirculationWorkspace()
    XP = GREBClimate.xghost
    cw_vec = (:T1h, :dTxh, :term_north, :term_south)
    cw_ghosted = (:T1h, :X_work, :wz_ghost)
    for f in fieldnames(CirculationWorkspace)
        v = getfield(cw, f)
        n = f in cw_ghosted ? XP : X
        @test eltype(v) === Float32
        @test size(v) == (f in cw_vec ? (n,) : (n, Y))
        @test all(iszero, v)
    end

    # MonthlyAccumulator: every field (xdim, ydim). There is no `count`
    # field - output! divides by cjday_mon[mon] * ndt_days.
    ma = MonthlyAccumulator()
    for f in fieldnames(MonthlyAccumulator)
        v = getfield(ma, f)
        @test size(v) == (X, Y) && eltype(v) === Float32 && all(iszero, v)
    end

end

@testset "derive_fields! follows the input maps" begin
    fields = synthetic_fields()
    quiet(() -> init_model!(resolve(preset(:full_model)), fields))
    @test fields.cap_surf[60, 10] == GREBClimate.cap_ocean * fields.mldclim[60, 10, 1]

    fields.z_topo[60, 10] = 500.0f0          # an ocean cell becomes land
    fields.mldclim[70, 20, :] .= 80.0f0
    fields.uclim[1, 1, 1] = -3.0f0
    fields.vclim[2, 2, 2] = 4.0f0
    GREBClimate.derive_fields!(fields, Processes())

    @test fields.cap_surf[60, 10] == GREBClimate.cap_land
    @test fields.wz_air[60, 10] == exp(-500.0f0 / GREBClimate.z_air)
    @test fields.wz_vapor[60, 10] == exp(-500.0f0 / GREBClimate.z_vapor)
    @test fields.z_ocean[70, 20] == 240.0f0
    @test fields.cap_surf[70, 20] == GREBClimate.cap_ocean * 80.0f0
    @test (fields.uclim_neg[1, 1, 1], fields.uclim_pos[1, 1, 1]) == (-3.0f0, 0.0f0)
    @test (fields.vclim_neg[2, 2, 2], fields.vclim_pos[2, 2, 2]) == (0.0f0, 4.0f0)
    @test fields.uclim_pos .+ fields.uclim_neg == fields.uclim
    @test fields.vclim_pos .+ fields.vclim_neg == fields.vclim
    @test all(>=(0), fields.uclim_pos) && all(<=(0), fields.uclim_neg)
end
