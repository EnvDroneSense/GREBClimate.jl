# Physics kernels: tendencies!, hydro!, SWradiation!, diffusion/advection/circulation.

@testset "tendencies! Q_sens honors the atmosphere switch" begin
    # Q_sens = ct_sens * (Ta - Ts), checked against a hand-computed value.
    fields = ClimateFields()
    state = ModelState()
    ws = ModelWorkspace()
    slice = 1
    r = resolve(preset(:full_model))

    Ts = fill(288.0, GREBClimate.xdim, GREBClimate.ydim)
    Ta = fill(280.0, GREBClimate.xdim, GREBClimate.ydim)
    To = fill(285.0, GREBClimate.xdim, GREBClimate.ydim)
    q = fill(0.006, GREBClimate.xdim, GREBClimate.ydim)

    tend = tendencies!(340.0, Ts, Ta, To, q, fields, state, ws, slice, r)
    @test tend.Q_sens ≈ GREBClimate.ct_sens .* (Ta .- Ts)
    @test all(isfinite, tend.SW)
    @test all(isfinite, tend.LW_surf)
    @test all(isfinite, tend.dTa_crcl)
    @test all(isfinite, tend.dq_crcl)

    r_off = resolve(preset(:full_model; processes = (atmosphere = false,)))
    tend_off = tendencies!(340.0, Ts, Ta, To, q, fields, state, ws, slice, r_off)
    @test all(iszero, tend_off.Q_sens)
end

@testset "land is topography above 0 m, in every kernel" begin
    @test GREBClimate.is_land.([-1.0f0, 0.0f0, 0.1f0]) == [false, false, true]

    # At 270 K land and ocean give different ice cover; a cell at exactly 0 m
    # goes with the ocean
    fields = ClimateFields()
    fields.z_topo[1, 1], fields.z_topo[2, 1], fields.z_topo[3, 1] = -1.0f0, 0.0f0, 1.0f0
    sw = SWradiation!(fill(270.0f0, X, Y), fields, ModelState(), 1, Processes(), ModelWorkspace())
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
        fields.u_clim_neg[i, k, it] = 0.5 + 0.0001 * i
        fields.u_clim_pos[i, k, it] = 0.3 + 0.0001 * k
        fields.v_clim_neg[i, k, it] = 0.4 + 0.0002 * i
        fields.v_clim_pos[i, k, it] = 0.2 + 0.0002 * k
    end
    ws = ModelWorkspace()
    slice = 1
    p = Processes()

    test_is = [1, 2, 3, 50, 94, 95, 96]
    test_ks = [1, 11, 48]

    diffusion!(T1, GREBClimate.z_air, fields, ws)
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
    diffusion!(Float64.(T1), GREBClimate.z_air, fields, ws)
    @test ws.dX_diff == dX_diff_f32
    diffusion!(view(T1, :, :), GREBClimate.z_air, fields, ws)
    @test ws.dX_diff == dX_diff_f32

    advection!(T1, GREBClimate.z_air, fields, ws, slice, p)
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
    circulation!(T1, GREBClimate.z_air, dX_out, fields, ws, slice, p)
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
    @test_throws ErrorException hydro!(Ts, q, ClimateFields(), 1, Processes(), h, ModelWorkspace())
end

@testset "hydro! fitted rain: dq_rain and Q_lat_air values, with no limit on rain" begin
    fields = constant_fields(z_topo = 1.0)
    init_model!(resolve(preset(:full_model)), fields)
    # c_q is large enough that any per-step limit on rain inside hydro! would
    # change the result: the fitted scheme applies none
    h = ResolvedHydrology(:fitted, :original, 1000, 0, 0, 0)

    Ts = fill(290.0f0, GREBClimate.xdim, GREBClimate.ydim)
    q = fill(0.008f0, GREBClimate.xdim, GREBClimate.ydim)
    slice = 1
    ws = ModelWorkspace()
    result = hydro!(Ts, q, fields, slice, Processes(), h, ws)

    expected_dq_rain = h.c_q * GREBClimate.cq_rain * q[1, 1]
    @test isapprox(result.dq_rain[1, 1], expected_dq_rain; rtol = 1e-5)
    @test isapprox(result.Q_lat_air[1, 1], -expected_dq_rain * GREBClimate.cq_latent * GREBClimate.r_qviwv; rtol = 1e-5)
end

@testset "hydro! evaporation: latent heat flux of the four schemes over land and ocean" begin
    G = GREBClimate
    Ts = fill(290.0f0, X, Y)
    q = fill(0.008f0, X, Y)
    qsat(T, wz) = 3.75e-3 * exp(17.08085 * (T - 273.15) / (T - 273.15 + 234.175)) * wz
    bulk = G.cq_latent * G.ρ_air * G.ce
    for land in (true, false)
        fields = constant_fields(z_topo = land ? 1.0 : -1.0, swet = 0.4, u = 3.0, v = 4.0, ws = 6.0)
        G.derive_fields!(fields, Processes())
        wz = fields.wz_air[1, 1]
        # :original and :original_gust use the wind components, the other
        # two the wind-speed climatology; :skin takes the saturation humidity
        # 5 K (land) or 1 K (ocean) above the surface temperature
        expected = (
            original = (0.008 - qsat(290, wz)) * sqrt(25 + (land ? 4 : 9)) * bulk * 0.4,
            original_gust = (0.008 - qsat(290, wz)) * sqrt(25 + (land ? 4 + 144 : 9 + 50.41)) *
                            bulk * (land ? 0.04 : 0.73) * 0.4,
            skin = (0.008 - qsat(290 + (land ? 5 : 1), wz)) * sqrt(36 + (land ? 132.25 : 29.16)) *
                   bulk * (land ? 0.25 : 0.58) * 0.4,
            skin_gust = (0.008 - qsat(290, wz)) * sqrt(36 + (land ? 81 : 16)) * bulk * (land ? 0.56 : 0.79) * 0.4,
        )
        for scheme in keys(expected)
            h = resolve(preset(:full_model; hydrology = (evaporation = scheme,))).hydrology
            result = hydro!(Ts, q, fields, 1, Processes(), h, ModelWorkspace())
            @test isapprox(result.Q_lat[1, 1], expected[scheme]; rtol = 1e-5)
            @test isapprox(result.dq_eva[1, 1], -expected[scheme] / G.cq_latent / G.r_qviwv; rtol = 1e-5)
        end
    end
end

@testset "SWradiation!: ice cover, albedo ramp and absorbed shortwave" begin
    G = GREBClimate
    fields = constant_fields(z_topo = -1.0)   # ocean, cloud cover 0.5
    fields.z_topo[1:3, 2] .= 1.0f0            # three land cells
    fields.glacier[4, 1] = 1.0f0              # a warm cell under a glacier
    fields.sw_solar .= 400.0f0
    # Fully ice-covered, the middle of the ramp, ice-free: ocean, then land
    Ts = fill(280.0f0, X, Y)
    Ts[1:2, 1] .= (260.0f0, (G.To_ice1 + G.To_ice2) / 2)
    Ts[1:2, 2] .= (250.0f0, (G.Tl_ice1 + G.Tl_ice2) / 2)
    state = ModelState()
    run(p) = map(copy, SWradiation!(Ts, fields, state, 1, p, ModelWorkspace()))

    a_atmos = 0.5 * G.a_cloud
    combined(a_surf) = a_surf + a_atmos - a_surf * a_atmos
    ramp = combined.(G.a_no_ice .+ G.da_ice .* [1, 0.5, 0])
    sw = run(Processes())
    for row in 1:2
        @test isapprox(sw.ice_cover[1:3, row], [1, 0.5, 0]; atol = 1e-4)
        @test isapprox(sw.albedo[1:3, row], ramp; rtol = 1e-4)
        @test isapprox(sw.SW[1:3, row], 400 .* (1 .- ramp); rtol = 1e-4)
    end
    # A glacier has the ice albedo at any temperature; it is not sea ice
    @test isapprox(sw.albedo[4, 1], ramp[1]; rtol = 1e-5) && sw.ice_cover[4, 1] == 0

    # The solar multiplier scales the absorbed flux
    state.sw_solar_forcing = 1.1f0
    @test isapprox(run(Processes()).SW, 1.1f0 .* sw.SW; rtol = 1e-5)
    state.sw_solar_forcing = 1.0f0

    # Without the ice-albedo feedback every surface is ice-free for the albedo
    off = run(Processes(ice_albedo = false))
    @test all(a -> isapprox(a, ramp[3]; rtol = 1e-5), off.albedo)
    @test off.ice_cover == sw.ice_cover

    # Exactly on a threshold: full ice at the lower one, none at the upper one
    for (T1, T2, inv) in ((G.To_ice1, G.To_ice2, G.inv_To_ice_range), (G.Tl_ice1, G.Tl_ice2, G.inv_Tl_ice_range))
        @test GREBClimate.@ice_ramp(T1, T1, T2, inv) === 1.0f0
        @test GREBClimate.@ice_ramp(T2, T1, T2, inv) === 0.0f0
        @test GREBClimate.@ice_ramp(prevfloat(T2), T1, T2, inv) > 0
    end
end

@testset "LWradiation!: emissivity and the longwave fluxes" begin
    G = GREBClimate
    fields = constant_fields(z_topo = 1500.0)   # cloud cover 0.5, Ts_clim 280 K
    G.derive_fields!(fields, Processes())
    Ts, Ta, q = fill(288.0f0, X, Y), fill(280.0f0, X, Y), fill(0.006f0, X, Y)
    run(co2; p = Processes(), f = fields) =
        map(copy, LWradiation!(Ts, Ta, q, co2, f, 1, p, ModelWorkspace()))
    lw = run(340.0f0)

    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = Float64.(G.emissivity_fit)
    wz = exp(-1500 / G.z_air)
    co2, vapor = wz * 340, wz * G.r_qviwv * 0.006
    clear = p4 * log(p1 * co2 + p2 * vapor + p3) + p7 + p5 * log(p1 * co2 + p3) + p6 * log(p2 * vapor + p3)
    em = (p8 - 0.5) / p9 * (clear - p10) + p10
    @test isapprox(lw.em[1, 1], em; rtol = 1e-4)
    @test isapprox(lw.LW_surf[1, 1], -G.σ * 288.0^4; rtol = 1e-5)
    T_rad = 280 + (-0.16 * 280 - 5)             # air temperature plus the radiation offset
    @test isapprox(lw.LW_down[1, 1], -em * G.σ * T_rad^4; rtol = 1e-4)
    @test lw.LW_up == lw.LW_down

    # More CO2 raises the emissivity; a cell with half the share of doubled
    # CO2 has the control emissivity
    @test run(680.0f0).em[1, 1] > lw.em[1, 1]
    half = deepcopy(fields)
    half.co2_part .= 0.5f0
    @test isapprox(run(680.0f0; f = half).em, lw.em; rtol = 1e-6)

    # Without an atmosphere there is no back radiation; LW_up keeps its value
    off = run(340.0f0; p = Processes(atmosphere = false))
    @test all(iszero, off.LW_down) && off.LW_up == lw.LW_up
end

@testset "seaice!: surface heat capacity along the ice ramp" begin
    G = GREBClimate
    fields = constant_fields(z_topo = -1.0)   # ocean, 50 m mixed layer
    fields.z_topo[4, 1] = 1.0f0               # land
    fields.glacier[5, 1] = 1.0f0              # ocean under a glacier
    G.derive_fields!(fields, Processes())
    open_ocean = G.cap_ocean * 50
    # Fully ice-covered, the middle of the ramp, then ice-free
    Ts = fill(280.0f0, X, Y)
    Ts[1:2, 1] .= (260.0f0, (G.To_ice1 + G.To_ice2) / 2)
    slice = 1

    seaice!(Ts, fields, slice, Processes())
    @test isapprox(fields.cap_surf[1:5, 1],
                   [G.cap_land, (G.cap_land + open_ocean) / 2, open_ocean, G.cap_land, G.cap_land]; rtol = 1e-4)

    # Without the ice-albedo feedback the ocean keeps its open-water capacity
    seaice!(Ts, fields, slice, Processes(ice_albedo = false))
    @test isapprox(fields.cap_surf[1:5, 1], [open_ocean, open_ocean, open_ocean, G.cap_land, G.cap_land]; rtol = 1e-6)

    # Without an ocean the kernel leaves the capacity alone
    fields.cap_surf .= 7.0f0
    seaice!(Ts, fields, slice, Processes(ocean = :none))
    @test all(==(7.0f0), fields.cap_surf)
end

@testset "deep_ocean!: entrainment, detrainment and turbulent mixing" begin
    G = GREBClimate
    fields = constant_fields(z_topo = -1.0)   # ocean, 50 m mixed layer at every step
    fields.z_topo[4, 1] = 1.0f0               # land
    fields.mld_clim[1, 1, N] = 40.0f0          # the step before step 1: the layer deepens by 10 m
    fields.mld_clim[2, 1, N] = 60.0f0          # ... or shoals by 10 m
    G.derive_fields!(fields, Processes())     # deep-ocean depth: 3 x the deepest mixed layer
    Ts = fill(290.0f0, X, Y)
    Ts[3, 1] = 260.0f0                        # under sea ice
    To = fill(280.0f0, X, Y)
    ws = ModelWorkspace()
    r = deep_ocean!(Ts, To, fields, 1, Processes(), ws)
    turb, mix = G.turb_coeff, G.c_effmix
    close(a, b) = isapprox(a, b; rtol = 1e-4)

    # Deepening: the surface layer takes up deep water; below, turbulent mixing only
    @test close(r.dT_ocean[1, 1], mix * (10 / 50) * (280 - 290) + turb * (280 - 290) / 50)
    @test close(r.dTo[1, 1], turb * (290 - 280) / (150 - 50))
    # Shoaling: the water left behind joins the deep ocean (180 m deep here)
    @test close(r.dTo[2, 1], mix * (10 / (180 - 50)) * (290 - 280) + turb * (290 - 280) / (180 - 50))
    @test close(r.dT_ocean[2, 1], turb * (280 - 290) / 50)
    # Under sea ice: mixing against the temperature at which the ice ends
    @test close(r.dT_ocean[3, 1], turb * (280 - G.To_ice2) / 50)
    @test close(r.dTo[3, 1], turb * (G.To_ice2 - 280) / (150 - 50))
    # Land
    @test r.dT_ocean[4, 1] == 0 && r.dTo[4, 1] == 0

    # Only the full ocean has a deep layer
    for ocean in (:mixed_layer, :none)
        off = deep_ocean!(Ts, To, fields, 1, Processes(ocean = ocean), ws)
        @test all(iszero, off.dT_ocean) && all(iszero, off.dTo)
    end
end

# The allocation budget for SWradiation! (and every other kernel) lives in
# test_invariants.jl, alongside the return-type checks.

@testset "hydrology = :none: humidity stays put, whatever the humidity correction" begin
    # Evaporation, rain and transport are zero with the water cycle off, so the
    # only term left to move q is the flux correction. SpinUp(0) keeps the one
    # set here.
    fields = synthetic_fields()
    fields.q_flux_correction .= 1.0f-4
    cfg = preset(:full_model; processes = (hydrology = :none,), corrections = SpinUp(0))
    q_ini = quiet(() -> init_model!(resolve(cfg), deepcopy(fields))).q_ini
    result = quiet() do
        greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = "", fields = fields, allow_uninitialized = true)
    end
    for rec in result.ctrl
        @test rec.q == q_ini
    end
end
