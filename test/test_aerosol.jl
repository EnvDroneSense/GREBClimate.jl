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
    state = ModelState()
    GREBClimate.apply_shortwave_addons!(state, cfg, 1950, 1)
    @test all(==(1.0f0), state.sw_solar_forcing)

    cfg.solar_scenario = Dict(1950 => 0.998f0)
    GREBClimate.apply_shortwave_addons!(state, cfg, 1950, 1)
    @test all(==(0.998f0), state.sw_solar_forcing)
    @test_throws ErrorException GREBClimate.apply_shortwave_addons!(state, cfg, 1951, 1)
    @test_throws ErrorException GREBClimate.apply_shortwave_addons!(state, cfg, 1950.5, 1)
end

@testset "apply_shortwave_addons! applies the aerosol and composes exactly with the solar table and solar experiments" begin
    cfg = create_experiment_config(:full_model)
    @test cfg.aerosol === nothing
    state = ModelState()
    cfg.aerosol = AerosolScenario()
    GREBClimate.apply_shortwave_addons!(state, cfg, 2000, 1)
    @test all(==(1.0f0), state.sw_solar_forcing)          # empty scenario: exact no-op

    cfg.aerosol = AerosolScenario(eruptions = [Eruption(2000.0, :tropical; peak_aod = 0.1)])
    # step 150 of the year is t = 2000 + 149/730
    GREBClimate.apply_shortwave_addons!(state, cfg, 2000, 150)
    only_aero = copy(state.sw_solar_forcing)
    @test all(<(1.0f0), only_aero)

    # Solar first, then aerosol; the product is exact.
    cfg.solar_scenario = Dict(2000 => 0.99f0)
    state.sw_solar_forcing .= 1.0f0
    GREBClimate.apply_shortwave_addons!(state, cfg, 2000, 150)
    @test state.sw_solar_forcing == 0.99f0 .* only_aero

    # Composes the same way with an experiment that already modulates the
    # multiplier: core/model.jl broadcasts forcing()'s scalar into the vector,
    # then the add-ons multiply into it.
    cyc = create_experiment_config(:solar_cycle_11yr)
    cyc.aerosol = cfg.aerosol
    scalar = forcing(150, 2000, cyc, ClimateFields(), zeros(Float64, X, Y, 1)).sw_solar_forcing
    state.sw_solar_forcing .= scalar
    GREBClimate.apply_shortwave_addons!(state, cyc, 2000, 150)
    @test state.sw_solar_forcing == scalar .* only_aero

    # The per-timestep path allocates nothing.
    GREBClimate.apply_shortwave_addons!(state, cfg, 2000, 150)
    @test (@allocated GREBClimate.apply_shortwave_addons!(state, cfg, 2000, 150)) == 0
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
    @test h(4.0) / h(3.0) ≈ exp(-1.0) rtol = 1e-3
end

@testset "aerosol chain: step response is 0 to 1 and its derivative is the impulse response over tau_decay" begin
    sc = AerosolScenario()
    tr, td = sc.tau_rise, sc.tau_decay
    @test GREBClimate.chain_step(0.0, tr, td) == 0.0
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

@testset "delta_eddington: limits, energy conservation, published closed form, monotonicity" begin
    de = GREBClimate.delta_eddington
    @test de(0.0, 1.0, 0.7, 0.5) == (0.0, 1.0)

    R, T = de(0.3, 0.0, 0.7, 0.5)
    @test R ≈ 0.0 atol = 1e-12
    @test T ≈ exp(-0.3 / 0.5)

    # Conservative scattering: energy is conserved, and R matches the
    # independent closed form of Meador and Weaver (1980) eq. 24 in the
    # delta-scaled quantities.
    for g in (0.5, 0.7, 0.85), mu0 in (0.3, 0.5, 0.8), tau in (0.01, 0.1, 0.5, 2.0)
        gs = g / (1 + g)
        g1 = (3 - 3gs) / 4
        g3 = (2 - 3gs * mu0) / 4
        ts = (1 - g^2) * tau
        R24 = (g1 * ts + (g3 - g1 * mu0) * (1 - exp(-ts / mu0))) / (1 + g1 * ts)
        R, T = de(tau, 1.0, g, mu0)
        @test R ≈ R24 rtol = 1e-6
        @test R + T ≈ 1.0 atol = 1e-6
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

@testset "aerosol_transmission!: exactly 1 before and long after, tropics dimmest for a tropical eruption, multiplies in place" begin
    sc = AerosolScenario(eruptions = [Eruption(2000.0, :tropical; peak_aod = 0.1)])
    buf = zeros(Float32, Y)
    mult = ones(Float32, Y)
    GREBClimate.aerosol_transmission!(mult, buf, sc, 1999.0)
    @test all(==(1.0f0), mult)
    GREBClimate.aerosol_transmission!(mult, buf, sc, 2050.0)   # dimming below Float32 resolution
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
end

@testset "tau_scale: calibrated by default, scales only the optical depth the transmission sees" begin
    e = Eruption(2000.0, :tropical; peak_aod = 0.1)
    @test AerosolScenario().tau_scale == GREBClimate.AEROSOL_TAU_SCALE
    buf = zeros(Float32, Y)
    trans(scale) = GREBClimate.aerosol_transmission!(ones(Float32, Y), buf,
        AerosolScenario(eruptions = [e], tau_scale = scale), 2000.42)

    for scale in (1.0, GREBClimate.AEROSOL_TAU_SCALE)
        mult = trans(scale)
        GREBClimate.aerosol_optical_depth!(buf, AerosolScenario(eruptions = [e]), 2000.42)
        @test mult == [Float32(GREBClimate.delta_eddington(scale * Float64(buf[j]), 1.0, 0.7, 0.5)[2])
                       for j in 1:Y]
    end
    @test all(==(1.0f0), trans(0.0))
    @test_throws ArgumentError AerosolScenario(tau_scale = -1.0)
    @test_throws ArgumentError AerosolScenario(tau_scale = NaN)
end

# Writes a forcing-series file of kind `aerosol_aod`; keywords override or `drop` entries.
function _write_aerosol(path; years, lat, aod, source = "test", version = 2, kind = "aerosol_aod",
                        calendar = "greb_365", units = "1", drop = ())
    entries = ("format_version" => version, "kind" => kind, "calendar" => calendar,
               "time" => years, "lat" => lat, "values" => aod, "units" => units, "source" => source)
    GREBClimate.jldopen(path, "w") do f
        for (k, v) in entries
            k in drop || (f[k] = v)
        end
    end
    return path
end

@testset "load_aerosol_series interpolates in time and latitude, and is zero outside the record" begin
    with_tempdir() do dir
        path = _write_aerosol(joinpath(dir, "aod.jld2"); years = [1991.0, 1992.0],
                              lat = [-60.0, 0.0, 60.0], aod = [0.0 0.2; 0.0 0.1; 0.0 0.2],
                              source = "unit test")
        # This record ends at 0.2, so loading it warns about the abrupt end.
        series = @test_logs (:warn, r"aod\.jld2.*1992.*0\.2") load_aerosol_series(path)
        @test series.source == "unit test"
        sc = AerosolScenario(series = series)

        aod = zeros(Float32, Y)
        GREBClimate.aerosol_optical_depth!(aod, sc, 1991.5)      # halfway in time
        @test (@allocated GREBClimate.aerosol_optical_depth!(aod, sc, 1991.5)) == 0
        for j in 1:Y
            phi = Float64(GREBClimate.lat_grid[j])
            band = phi <= -60 ? 0.2 : phi >= 60 ? 0.2 :
                   phi <= 0 ? 0.2 + (0.1 - 0.2) * (phi + 60) / 60 :
                              0.1 + (0.2 - 0.1) * phi / 60
            @test aod[j] ≈ 0.5 * band rtol = 1e-5
        end

        GREBClimate.aerosol_optical_depth!(aod, sc, 1990.0)
        @test all(iszero, aod)
        GREBClimate.aerosol_optical_depth!(aod, sc, 1993.0)
        @test all(iszero, aod)

        # A series and an eruption over the same period add.
        e = Eruption(1991.5, :tropical; peak_aod = 0.1)
        both = AerosolScenario(series = series, eruptions = [e])
        only_e = AerosolScenario(eruptions = [e])
        a_series, a_e, a_both = (zeros(Float32, Y) for _ in 1:3)
        GREBClimate.aerosol_optical_depth!(a_series, sc, 1991.7)
        GREBClimate.aerosol_optical_depth!(a_e, only_e, 1991.7)
        GREBClimate.aerosol_optical_depth!(a_both, both, 1991.7)
        @test a_both ≈ a_series .+ a_e

        # Exact data years return their column; just past the last year the record is zero.
        aod0 = zeros(Float32, Y)
        GREBClimate.aerosol_optical_depth!(aod0, sc, 1991.0)
        @test all(iszero, aod0)
        aodN = zeros(Float32, Y)
        GREBClimate.aerosol_optical_depth!(aodN, sc, 1992.0)
        @test aodN == series.aod[:, 2]
        GREBClimate.aerosol_optical_depth!(aodN, sc, 1992.0 + 1e-9)
        @test all(iszero, aodN)

        # A record that ends below the threshold loads silently and unaltered.
        fades = _write_aerosol(joinpath(dir, "fades.jld2"); years = [1991.0, 1992.0, 1993.0],
                               lat = [-60.0, 0.0, 60.0],
                               aod = [0.0 0.2 0.0; 0.0 0.1 0.0; 0.0 0.2 0.001])
        faded = @test_logs load_aerosol_series(fades)
        @test faded.aod[1, end] == 0.0f0

        good = (years = [1991.0, 1992.0], lat = [-60.0, 0.0, 60.0], aod = zeros(3, 2))
        for (name, kw) in (
                ("nokey.jld2",     (; good..., drop = ("values",))),
                ("version.jld2",   (; good..., version = 1)),
                ("noversion.jld2", (; good..., drop = ("format_version",))),
                ("kind.jld2",      (; good..., kind = "tsi")),
                ("calendar.jld2",  (; good..., calendar = "gregorian")),
                ("units.jld2",     (; good..., units = "W/m2")),
                ("nolat.jld2",     (; good..., drop = ("lat",), aod = zeros(2))),
                ("shape.jld2",     (; good..., aod = zeros(3, 3))),
                ("unsorted.jld2",  (; good..., years = [1992.0, 1991.0])),
                ("latorder.jld2",  (; good..., lat = [60.0, 0.0, -60.0])),
                ("nan.jld2",       (; good..., aod = [0.0 0.0; NaN 0.0; 0.0 0.0])),
                ("empty.jld2",     (; years = Float64[], lat = [-60.0, 0.0, 60.0], aod = zeros(3, 0))))
            @test_throws ErrorException load_aerosol_series(_write_aerosol(joinpath(dir, name); kw...))
        end
        @test_throws ErrorException load_aerosol_series(joinpath(dir, "missing.jld2"))
    end
end
