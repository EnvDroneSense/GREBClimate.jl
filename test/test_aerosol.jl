# Shortwave add-ons: the solar series and the stratospheric aerosol kernel.

@testset "load_solar_series returns multipliers relative to the table mean or a reference" begin
    with_tempdir() do dir
        path = joinpath(dir, "tsi.txt")
        write(path, "# year TSI\n1990 1360.0\n\n1991 1362.0\n1992 1361.0\n")
        s = load_solar_series(path)
        @test s[1991] ≈ Float32(1362 / 1361)
        @test s[1990] ≈ Float32(1360 / 1361)
        s2 = load_solar_series(path; reference = 1365.0)
        @test s2[1991] ≈ Float32(1362 / 1365)

        bad = joinpath(dir, "bad.txt")
        write(bad, "1990\n")
        @test_throws ErrorException load_solar_series(bad)
        @test_throws ErrorException load_solar_series(joinpath(dir, "missing.txt"))
        empty_path = joinpath(dir, "empty.txt")
        write(empty_path, "# nothing\n")
        @test_throws ErrorException load_solar_series(empty_path)
        @test_throws ErrorException load_solar_series(path; reference = -1.0)

        gap = joinpath(dir, "gap.txt")
        write(gap, "1990 1360.0\n1992 1361.0\n")
        @test_throws "gap" load_solar_series(gap)
        dup = joinpath(dir, "dup.txt")
        write(dup, "1990 1360.0\n1990 1361.0\n")
        @test_throws "more than once" load_solar_series(dup)
    end
end

@testset "apply_shortwave_addons! is an exact no-op unconfigured and applies the solar table otherwise" begin
    cfg = create_experiment_config(:full_model)
    mult = ones(Float32, Y)
    GREBClimate.apply_shortwave_addons!(mult, cfg, 1950, 1)
    @test all(==(1.0f0), mult)

    cfg.solar_scenario = Dict(1950 => 0.998f0)
    GREBClimate.apply_shortwave_addons!(mult, cfg, 1950, 1)
    @test all(==(0.998f0), mult)
    @test_throws ErrorException GREBClimate.apply_shortwave_addons!(mult, cfg, 1951, 1)
    @test_throws ErrorException GREBClimate.apply_shortwave_addons!(mult, cfg, 1950.5, 1)
end

const _W = cosd.(Float64.(GREBClimate.lat_grid))
_gm(a) = sum(a .* _W) / sum(_W)

@testset "aerosol chain: impulse peaks at the analytic peak time and decays with tau_decay" begin
    sc = AerosolScenario()
    tr, td = sc.tau_rise, sc.tau_decay
    h(t) = GREBClimate.chain_impulse(t, tr, td)
    tpk = GREBClimate.chain_peak_time(tr, td)
    # Crowley and Unterman (2013): the tropical Pinatubo peak came about 5 months after the eruption.
    @test 12 * tpk ≈ 5.0 atol = 0.1
    @test maximum(h.(range(0.0, 3.0; length = 3001))) ≈ h(tpk) rtol = 1e-6
    @test h(0.0) == 0.0
    @test h(-0.5) == 0.0
    @test h(4.0) / h(3.0) ≈ exp(-1.0) rtol = 1e-3
end

@testset "aerosol chain: step response is 0 to 1 and its derivative is the impulse response over tau_decay" begin
    sc = AerosolScenario()
    tr, td = sc.tau_rise, sc.tau_decay
    @test GREBClimate.chain_step(0.0, tr, td) == 0.0
    @test GREBClimate.chain_step(-1.0, tr, td) == 0.0
    @test GREBClimate.chain_step(60.0, tr, td) ≈ 1.0
    t, d = 0.3, 1e-5
    slope = (GREBClimate.chain_step(t + d, tr, td) - GREBClimate.chain_step(t - d, tr, td)) / 2d
    @test slope ≈ GREBClimate.chain_impulse(t, tr, td) / td rtol = 1e-4
end

@testset "Eruption and SustainedInjection validate their arguments" begin
    @test Eruption(1991.45, :tropical; tg_s = 10.0).peak_aod ≈ 10.0 * GREBClimate.AOD_PER_TG_S
    @test Eruption(1991.45, :tropical; peak_aod = 0.12).peak_aod == 0.12
    @test_throws ArgumentError Eruption(1991.45, :tropical)
    @test_throws ArgumentError Eruption(1991.45, :tropical; tg_s = 10.0, peak_aod = 0.1)
    @test_throws ArgumentError Eruption(1991.45, :polar; peak_aod = 0.1)
    @test_throws ArgumentError Eruption(1991.45, :tropical; peak_aod = -0.1)
    @test_throws ArgumentError SustainedInjection(2020.0, 2010.0, :tropical, 0.05)
    @test_throws ArgumentError SustainedInjection(2010.0, 2020.0, :tropical, -0.05)
    @test SustainedInjection(2010.0, Inf, :tropical, 0.05).stop_year == Inf
    @test_throws ArgumentError AerosolScenario(tau_rise = 2.0, tau_decay = 1.0)
    @test_throws ArgumentError AerosolScenario(ssa = 1.5)
end

@testset "class profiles are area-normalised, symmetric where expected, and match the sourced tropical shape" begin
    for cls in GREBClimate.AEROSOL_CLASSES
        @test _gm(GREBClimate._CLASS_PROFILES[cls]) ≈ 1.0
    end
    trop = GREBClimate._CLASS_PROFILES[:tropical]
    @test trop ≈ reverse(trop)
    @test trop[Y ÷ 2] == maximum(trop)
    @test issorted(trop[1:Y ÷ 2])
    # Crowley and Unterman (2013): high-latitude values are about 80% of the tropical ones.
    @test trop[1] / trop[Y ÷ 2] ≈ 0.8
    nh = GREBClimate._CLASS_PROFILES[:nh_extratropical]
    sh = GREBClimate._CLASS_PROFILES[:sh_extratropical]
    @test nh ≈ reverse(sh)
    @test sum(nh[Y ÷ 2 + 1:end]) > sum(nh[1:Y ÷ 2])
end

@testset "aerosol_optical_depth!: global-mean peak, superposition, zero before, sustained injection" begin
    sc0 = AerosolScenario()
    tr, td = sc0.tau_rise, sc0.tau_decay
    for cls in GREBClimate.AEROSOL_CLASSES
        sc = AerosolScenario(eruptions = [Eruption(2000.0, cls; peak_aod = 0.1)])
        aod = zeros(Float32, Y)
        GREBClimate.aerosol_optical_depth!(aod, sc, 2000.0 + GREBClimate.chain_peak_time(tr, td))
        @test _gm(aod) ≈ 0.1 rtol = 1e-5
        GREBClimate.aerosol_optical_depth!(aod, sc, 1999.0)
        @test all(iszero, aod)
    end

    e1 = Eruption(2000.0, :tropical; peak_aod = 0.1)
    e2 = Eruption(2001.0, :nh_extratropical; peak_aod = 0.05)
    a1, a2, a12 = (zeros(Float32, Y) for _ in 1:3)
    GREBClimate.aerosol_optical_depth!(a1, AerosolScenario(eruptions = [e1]), 2001.7)
    GREBClimate.aerosol_optical_depth!(a2, AerosolScenario(eruptions = [e2]), 2001.7)
    GREBClimate.aerosol_optical_depth!(a12, AerosolScenario(eruptions = [e1, e2]), 2001.7)
    @test a12 ≈ a1 .+ a2

    sc = AerosolScenario(injections = [SustainedInjection(2010.0, 2020.0, :tropical, 0.05)])
    aod = zeros(Float32, Y)
    GREBClimate.aerosol_optical_depth!(aod, sc, 2015.0)
    @test _gm(aod) ≈ 0.05 * GREBClimate.chain_step(5.0, tr, td) rtol = 1e-5
    GREBClimate.aerosol_optical_depth!(aod, sc, 2025.0)
    @test _gm(aod) ≈ 0.05 * (GREBClimate.chain_step(15.0, tr, td) - GREBClimate.chain_step(5.0, tr, td)) rtol = 1e-3
    # termination: the anomaly decays with tau_decay once the source stops
    a2y, a3y = zeros(Float32, Y), zeros(Float32, Y)
    GREBClimate.aerosol_optical_depth!(a2y, sc, 2022.0)
    GREBClimate.aerosol_optical_depth!(a3y, sc, 2023.0)
    @test _gm(a3y) / _gm(a2y) ≈ exp(-1.0) rtol = 1e-2
end

@testset "delta_eddington: limits, energy conservation, thin-layer reflectance, published closed forms, monotonicity" begin
    de = GREBClimate.delta_eddington
    @test de(0.0, 1.0, 0.7, 0.5) == (0.0, 1.0)

    R, T = de(0.3, 0.0, 0.7, 0.5)
    @test R ≈ 0.0 atol = 1e-12
    @test T ≈ exp(-0.3 / 0.5)

    for tau in (0.05, 0.5, 2.0)
        R, T = de(tau, 1.0, 0.7, 0.5)
        @test R + T ≈ 1.0 atol = 1e-6
    end

    g, mu0 = 0.7, 0.5
    gs = g / (1 + g)
    gamma3 = (2 - 3gs * mu0) / 4
    R, _ = de(1e-4, 1.0, g, mu0)
    @test R / 1e-4 ≈ (1 - g^2) * gamma3 / mu0 rtol = 1e-3

    # Independent closed form for conservative scattering, Meador and Weaver
    # (1980) eq. 24, in the delta-scaled quantities.
    for g in (0.5, 0.7, 0.85), mu0 in (0.3, 0.5, 0.8), tau in (0.01, 0.1, 0.5, 2.0)
        gs = g / (1 + g)
        g1 = (3 - 3gs) / 4
        g3 = (2 - 3gs * mu0) / 4
        ts = (1 - g^2) * tau
        R24 = (g1 * ts + (g3 - g1 * mu0) * (1 - exp(-ts / mu0))) / (1 + g1 * ts)
        @test de(tau, 1.0, g, mu0)[1] ≈ R24 rtol = 1e-6
    end

    # Joseph et al. (1976) similarity relation (eq. 17b) for the delta scaling.
    for w in 0.1:0.2:0.9, g in 0.1:0.2:0.9
        f = g^2
        wp = (1 - f) * w / (1 - w * f)
        gp = g / (1 + g)
        @test (1 - wp) / (1 - wp * gp) ≈ (1 - w) / (1 - w * g)
    end

    taus = 0.01:0.01:1.0
    RT = [de(t, 1.0, 0.7, 0.5) for t in taus]
    @test issorted(first.(RT))
    @test issorted(last.(RT); rev = true)

    @test_throws DomainError de(200.0, 1.0, 0.7, 0.5)      # tau / mu0 = 400
    @test isfinite(de(100.0, 1.0, 0.7, 0.5)[2])            # tau / mu0 = 200, still safe
end

@testset "delta_eddington stays finite across ssa and at the k*mu0 = 1 singularity" begin
    de = GREBClimate.delta_eddington
    for ssa in 0.0:0.01:1.0
        R, T = de(0.2, ssa, 0.7, 1.0)
        @test isfinite(R) && isfinite(T)
        @test 0.0 <= T <= 1.0 + 1e-9
    end
    ssa, g = 0.5, 0.7
    f = g^2
    om = (1 - f) * ssa / (1 - ssa * f)
    gs = g / (1 + g)
    g1 = (7 - om * (4 + 3gs)) / 4
    g2 = -(1 - om * (4 - 3gs)) / 4
    mu_sing = 1 / sqrt(g1^2 - g2^2)
    R, T = de(0.2, ssa, g, mu_sing)
    @test isfinite(R) && isfinite(T)
    # The nudge must not be visible: the result agrees with values a small
    # distance either side of the singularity (the singularity is removable).
    for eps in (1e-3, -1e-3)
        Re, Te = de(0.2, ssa, g, mu_sing * (1 + eps))
        @test abs(Re - R) < 1e-3
        @test abs(Te - T) < 1e-3
    end
end

@testset "aerosol_transmission!: exactly 1 with no aerosol, tropics dimmest for a tropical eruption, multiplies in place" begin
    sc = AerosolScenario(eruptions = [Eruption(2000.0, :tropical; peak_aod = 0.1)])
    buf = zeros(Float32, Y)
    mult = ones(Float32, Y)
    GREBClimate.aerosol_transmission!(mult, buf, sc, 1999.0)
    @test all(==(1.0f0), mult)
    GREBClimate.aerosol_transmission!(mult, buf, sc, 2000.42)
    @test all(<(1.0f0), mult)
    @test mult[Y ÷ 2] == minimum(mult)

    half = fill(0.5f0, Y)
    GREBClimate.aerosol_transmission!(half, buf, sc, 1999.0)
    @test all(==(0.5f0), half)

    nh = AerosolScenario(eruptions = [Eruption(2000.0, :nh_extratropical; peak_aod = 0.1)])
    m = ones(Float32, Y)
    GREBClimate.aerosol_transmission!(m, buf, nh, 2000.42)
    @test sum(m[Y ÷ 2 + 1:end]) < sum(m[1:Y ÷ 2])

    # The per-timestep path allocates nothing.
    GREBClimate.aerosol_transmission!(m, buf, nh, 2000.42)
    @test (@allocated GREBClimate.aerosol_transmission!(m, buf, nh, 2000.42)) == 0
end
