# Physics kernels: tendencies!, hydro!, SWradiation!, diffusion/advection/circulation.

@testset "tendencies! Q_sens honors the atmosphere switch" begin
    # Q_sens = ct_sens * (Ta - Ts) is checkable directly against a
    # hand-computed value without reimplementing the rest of the
    # physics pipeline.
    fields = ClimateFields()
    state = ModelState()
    ws = CirculationWorkspace()
    ts = TimeState(1, 1)
    r = resolve(preset(:full_model))

    Ts = fill(288.0, GREBClimate.xdim, GREBClimate.ydim)
    Ta = fill(280.0, GREBClimate.xdim, GREBClimate.ydim)
    To = fill(285.0, GREBClimate.xdim, GREBClimate.ydim)
    q = fill(0.006, GREBClimate.xdim, GREBClimate.ydim)

    tend = tendencies!(340.0, Ts, Ta, To, q, fields, state, ws, ts, r)
    @test tend.Q_sens ≈ GREBClimate.ct_sens .* (Ta .- Ts)
    @test all(isfinite, tend.SW)
    @test all(isfinite, tend.LW_surf)
    @test all(isfinite, tend.dTa_crcl)
    @test all(isfinite, tend.dq_crcl)

    r_off = resolve(preset(:full_model; processes = (atmosphere = false,)))
    tend_off = tendencies!(340.0, Ts, Ta, To, q, fields, state, ws, ts, r_off)
    @test all(iszero, tend_off.Q_sens)
end

@testset "land is topography above 0 m, in every kernel" begin
    @test GREBClimate.is_land.([-1.0f0, 0.0f0, 0.1f0]) == [false, false, true]

    # At 270 K land and ocean give different ice cover; a cell at exactly 0 m
    # goes with the ocean
    fields = ClimateFields()
    fields.z_topo[1, 1], fields.z_topo[2, 1], fields.z_topo[3, 1] = -1.0f0, 0.0f0, 1.0f0
    sw = SWradiation!(fill(270.0f0, X, Y), fields, ModelState(), TimeState(1, 1), Processes(), CirculationWorkspace())
    @test sw.ice_cover[2, 1] == sw.ice_cover[1, 1]
    @test sw.ice_cover[3, 1] > sw.ice_cover[1, 1]
    @test sw.albedo[2, 1] == sw.albedo[1, 1]
end

@testset "diffusion!/advection!/circulation! per-cell snapshot (incl. date-line wraparound)" begin
    fields = ClimateFields()
    xdim_, ydim_ = GREBClimate.xdim, GREBClimate.ydim
    T1 = Float32[100.0 * i + k for i in 1:xdim_, k in 1:ydim_]
    wz = [1.0 + 0.001 * i - 0.0005 * k for i in 1:xdim_, k in 1:ydim_]
    fields.wz_air .= wz
    fields.wz_vapor .= wz
    for it in 1:GREBClimate.nstep_yr, k in 1:ydim_, i in 1:xdim_
        fields.uclim_neg[i, k, it] = 0.5 + 0.0001 * i
        fields.uclim_pos[i, k, it] = 0.3 + 0.0001 * k
        fields.vclim_neg[i, k, it] = 0.4 + 0.0002 * i
        fields.vclim_pos[i, k, it] = 0.2 + 0.0002 * k
    end
    ws = CirculationWorkspace()
    ts = TimeState(1, 1)
    p = Processes()

    test_is = [1, 2, 3, 50, 94, 95, 96]
    test_ks = [1, 11, 48]

    diffusion!(T1, GREBClimate.z_air, fields, ws, ts)
    dX_diff_ref = Dict(
        (1,1)=>4517.8355192140425, (2,1)=>3542.95440832927, (3,1)=>2652.6121358299374,
        (50,1)=>1.2269892658145531, (94,1)=>-2709.198844945152, (95,1)=>-3576.970530673063,
        (96,1)=>-4530.555031256096,
        (1,11)=>64.24343226619055, (2,11)=>32.16535999519279, (3,11)=>10.737877625759896,
        (50,11)=>0.0032154796768526627, (94,11)=>-10.710851102068029, (95,11)=>-32.17954307109926,
        (96,11)=>-64.44305854028968,
        (1,48)=>4410.353233774593, (2,48)=>3445.760591202237, (3,48)=>2567.190165268745,
        (50,48)=>1.1818789937009289, (94,48)=>-2623.7462653181888, (95,48)=>-3479.421880797481,
        (96,48)=>-4422.060857683985,
    )
    for k in test_ks, i in test_is
        @test isapprox(ws.dX_diff[i, k], dX_diff_ref[(i, k)]; atol=1e-3, rtol=1e-4)
    end

    # The public API accepts Float64 arrays and views. Those miss `to_ghosted!`'s
    # `Matrix{Float32}` memcpy fast path and take the generic fallback, which
    # the snapshots above never reach. These values are exactly representable in
    # Float32, so both paths must agree bit for bit.
    dX_diff_f32 = copy(ws.dX_diff)
    diffusion!(Float64.(T1), GREBClimate.z_air, fields, ws, ts)
    @test ws.dX_diff == dX_diff_f32
    diffusion!(view(T1, :, :), GREBClimate.z_air, fields, ws, ts)
    @test ws.dX_diff == dX_diff_f32

    advection!(T1, GREBClimate.z_air, fields, ws, ts, p)
    dX_adv_ref = Dict(
        (1,1)=>99.99281072836801, (2,1)=>37.62423619149657, (3,1)=>6.428217314852662,
        (50,1)=>-4.182755344492395, (94,1)=>11.771946283998448, (95,1)=>60.257791236566284,
        (96,1)=>157.2916794570926,
        (1,11)=>6.863756666482768, (2,11)=>3.295185364395947, (3,11)=>-0.2733860065994545,
        (50,11)=>-0.28795610784365167, (94,11)=>-0.30173415760057604, (95,11)=>5.231130705764574,
        (96,11)=>10.766167503970431,
        (1,48)=>99.42365199872042, (2,48)=>37.439173701522, (3,48)=>6.435048732922689,
        (50,48)=>-4.112511103625333, (94,48)=>11.462671629384038, (95,48)=>58.8110669109033,
        (96,48)=>153.56947120915834,
    )
    for k in test_ks, i in test_is
        @test isapprox(ws.dX_adv[i, k], dX_adv_ref[(i, k)]; atol=1e-3, rtol=1e-4)
    end

    dX_out = zeros(xdim_, ydim_)
    circulation!(T1, GREBClimate.z_air, dX_out, fields, ws, ts, p)
    dX_out_ref = Dict(
        (1,1)=>4824.155681112328, (2,1)=>4609.661409972268, (3,1)=>4395.402855867866,
        (50,1)=>-98.5576039612888, (94,1)=>-4172.300403641432, (95,1)=>-4369.7800972827745,
        (96,1)=>-4569.32742911384,
        (1,11)=>1447.9977087179682, (2,11)=>765.1034286356905, (3,11)=>274.82096986927763,
        (50,11)=>-6.684077032670757, (94,11)=>-252.25457239155912, (95,11)=>-576.3346453253889,
        (96,11)=>-1090.2114380116673,
        (1,48)=>4827.018868288203, (2,48)=>4605.028237288625, (3,48)=>4383.415249532827,
        (50,48)=>-94.78570775574462, (94,48)=>-4150.449604784975, (95,48)=>-4353.86760064195,
        (96,48)=>-4559.610385060042,
    )
    for k in test_ks, i in test_is
        @test isapprox(dX_out[i, k], dX_out_ref[(i, k)]; atol=1e-3, rtol=1e-4)
    end
end

@testset "hydro! errors on an unknown evaporation scheme" begin
    Ts = fill(290.0, GREBClimate.xdim, GREBClimate.ydim)
    q = fill(0.005, GREBClimate.xdim, GREBClimate.ydim)
    h = ResolvedHydrology(:fitted, :bogus, 1, 0, 0, 0)
    @test_throws ErrorException hydro!(Ts, q, ClimateFields(), TimeState(1, 1), Processes(), h, CirculationWorkspace())
end

@testset "hydro! :original_gust evaporation includes Fortran's carried-over +2.0²/+3.0² base term" begin
    mkfields(topo) = begin
        fields = ClimateFields()
        fields.z_topo .= topo
        fields.mldclim .= 50.0
        fields.Tclim .= 280.0
        fields.Toclim .= 285.0
        fields.qclim .= 0.006
        fields.cldclim .= 0.5
        fields.swetclim .= 1.0
        fields.uclim .= 0.0
        fields.vclim .= 0.0
        fields.omegaclim .= 0.0
        fields.omegastdclim .= 0.0
        fields.wsclim .= 0.0
        fields
    end
    r = resolve(preset(:full_model; hydrology = (evaporation = :original_gust,)))
    Ts = fill(290.0f0, GREBClimate.xdim, GREBClimate.ydim)
    q = fill(0.008f0, GREBClimate.xdim, GREBClimate.ydim)
    ts = TimeState(1, 1)

    for (topo, gust, coeff) in ((1.0, 4.0 + 144.0, 0.04), (-1.0, 9.0 + 50.41, 0.73))
        fields = mkfields(topo)
        init_model!(r, fields)
        ws = CirculationWorkspace()
        result = hydro!(Ts, q, fields, ts, Processes(), r.hydrology, ws)

        qs = 3.75e-3 * exp(17.08085 * (290.0 - 273.15) / (290.0 - 273.15 + 234.175)) * fields.wz_air[1, 1]
        expected = (q[1, 1] - qs) * sqrt(gust) * GREBClimate.cq_latent * GREBClimate.ρ_air * coeff * GREBClimate.ce * 1.0
        @test isapprox(result.Q_lat[1, 1], expected; rtol = 1e-5)
    end
end

@testset "hydro! doesn't apply an extra -0.9q clamp to dq_rain" begin
    fields = ClimateFields()
    fields.z_topo .= 1.0
    fields.mldclim .= 50.0
    fields.Tclim .= 280.0
    fields.Toclim .= 285.0
    fields.qclim .= 0.006
    fields.cldclim .= 0.5
    fields.swetclim .= 1.0
    fields.uclim .= 0.0
    fields.vclim .= 0.0
    fields.omegaclim .= 0.0
    fields.omegastdclim .= 0.0
    fields.wsclim .= 0.0
    init_model!(resolve(preset(:full_model)), fields)
    # The fitted scheme has no rain-limit clamp; a huge c_q makes the old -0.9q clamp fire
    h = ResolvedHydrology(:fitted, :original, 1000, 0, 0, 0)

    Ts = fill(290.0f0, GREBClimate.xdim, GREBClimate.ydim)
    q = fill(0.008f0, GREBClimate.xdim, GREBClimate.ydim)
    ts = TimeState(1, 1)
    ws = CirculationWorkspace()
    result = hydro!(Ts, q, fields, ts, Processes(), h, ws)

    expected_dq_rain = h.c_q * GREBClimate.cq_rain * q[1, 1]
    min_dq_that_would_have_clamped = -0.9 * q[1, 1] / GREBClimate.Δt
    @test expected_dq_rain < min_dq_that_would_have_clamped  # sanity: the old clamp would have fired
    @test isapprox(result.dq_rain[1, 1], expected_dq_rain; rtol = 1e-5)
    @test isapprox(result.Q_lat_air[1, 1], -expected_dq_rain * GREBClimate.cq_latent * GREBClimate.r_qviwv; rtol = 1e-5)
end

# The allocation budget for SWradiation! (and every other kernel) lives in
# test_invariants.jl, alongside the return-type checks.

@testset "hydrology = :none: humidity stays put, whatever the humidity correction" begin
    # Evaporation, rain and transport are zero with the water cycle off, so the
    # only term left to move q is the flux correction. SpinUp(0) keeps the one
    # set here.
    fields = synthetic_fields()
    fields.qF_correct .= 1.0f-4
    cfg = preset(:full_model; processes = (hydrology = :none,), corrections = SpinUp(0))
    q_ini = quiet(() -> init_model!(resolve(cfg), deepcopy(fields))).q_ini
    result = quiet() do
        greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = "", fields = fields, allow_uninitialized = true)
    end
    for rec in result.ctrl
        @test rec.q == q_ini
    end
end
