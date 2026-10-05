# Golden regression against a saved snapshot of a real-dataset run.

include(joinpath("data", "golden_reference.jl"))

@testset "golden regression: real dataset control+scenario run matches snapshot" begin
    # A 1-year control and 1-year scenario on the NCEP dataset, stored as
    # monthly global-mean Ts/Ta/q, plus a control year after a one-year spin-up.
    # The tolerances (1e-3 K, 1e-6 kg/kg) are about 30 times the drift measured
    # between runs; exact equality is checked by tools/validation/bit_identity.jl.
    # RUN_GOLDEN=0 skips it locally. CI has no dataset and always skips it, so
    # a golden break is local-red and CI-green.
    if !isdir(DATA_DIR)
        @test_skip "greb_input_data/ not present"
    elseif get(ENV, "RUN_GOLDEN", "1") == "0"
        @test_skip "RUN_GOLDEN=0"
    else
        fields = load_climatology(DATA_DIR; dataset = :ncep)
        result = quiet() do
            greb_model!(RunSpec(), preset(:full_model; corrections = Stored()); jld2_dir = DATA_DIR, fields = fields)
        end

        gmean(x) = sum(x) / length(x)
        summarize(rec) = (Ts = gmean(rec.Ts), Ta = gmean(rec.Ta), q = gmean(rec.q))

        @test length(result.ctrl) == length(GOLDEN_CTRL)
        @test length(result.scnr) == length(GOLDEN_SCNR)
        for (rec, ref) in zip(result.ctrl, GOLDEN_CTRL)
            s = summarize(rec)
            @test isapprox(s.Ts, ref.Ts; atol = 1e-3)
            @test isapprox(s.Ta, ref.Ta; atol = 1e-3)
            @test isapprox(s.q, ref.q; atol = 1e-6)
        end
        for (rec, ref) in zip(result.scnr, GOLDEN_SCNR)
            s = summarize(rec)
            @test isapprox(s.Ts, ref.Ts; atol = 1e-3)
            @test isapprox(s.Ta, ref.Ta; atol = 1e-3)
            @test isapprox(s.q, ref.q; atol = 1e-6)
        end

        # One spin-up year, then the control: exercises qflux_correction!.
        # Reuses `fields`, which greb_model! restores after the run above.
        flux_result = quiet() do
            greb_model!(RunSpec(ctrl = 1, scnr = 0), preset(:full_model; corrections = SpinUp(1));
                        jld2_dir = DATA_DIR, fields = fields)
        end
        @test length(flux_result.ctrl) == length(GOLDEN_FLUX)
        for (rec, ref) in zip(flux_result.ctrl, GOLDEN_FLUX)
            s = summarize(rec)
            @test isapprox(s.Ts, ref.Ts; atol = 1e-3)
            @test isapprox(s.Ta, ref.Ta; atol = 1e-3)
            @test isapprox(s.q, ref.q; atol = 1e-6)
        end
    end
end

# Area-weighted global mean of a monthly-record field, averaged over `recs`
function area_mean(recs, var = :Ts)
    w = cosd.(range(-88.125, 88.125; length = Y))
    return sum(sum(getfield(r, var) .* w') for r in recs) / (length(recs) * X * sum(w))
end

@testset "MSCM configuration reproduces the MSCM 2xCO2 response" begin
    # MSCM (Monash Simple Climate Model) database, 2xCO2 with every process on:
    # year-1 global-mean surface temperature response 0.594636 K.
    if !isdir(DATA_DIR)
        @test_skip "greb_input_data/ not present"
    else
        cfg = preset(:co2_double; processes = (moisture_convergence = false,), hydrology = mscm_hydrology())
        result = quiet() do
            greb_model!(RunSpec(ctrl = 1, scnr = 1), cfg; jld2_dir = DATA_DIR,
                        fields = load_climatology(DATA_DIR; dataset = :ncep))
        end
        @test isapprox(area_mean(result.scnr), MSCM_YEAR1; atol = 1e-3)
    end
end

@testset "mean-climate deconstruction: a switched-off process changes the control climate" begin
    # On computed corrections every configuration is pulled back to the
    # observed climate; on the stored ones the switch shows.
    if !isdir(DATA_DIR)
        @test_skip "greb_input_data/ not present"
    else
        fields = load_climatology(DATA_DIR; dataset = :ncep)
        control(processes) = quiet() do
            greb_model!(RunSpec(ctrl = 1, scnr = 0), preset(:decon_mean_climate; processes);
                        jld2_dir = DATA_DIR, fields)
        end.ctrl
        @test abs(area_mean(control((ocean = :none,))) - area_mean(control((;)))) > 0.1
    end
end

@testset "2xCO2 deconstruction: uniform humidity gives a finite response" begin
    # 0.0052 kg/kg is above saturation over cold, high ground. The fitted rain
    # scheme rains that out within the spin-up and the run goes non-finite;
    # the preset's MSCM physics does not.
    if !isdir(DATA_DIR)
        @test_skip "greb_input_data/ not present"
    else
        result = quiet() do
            greb_model!(RunSpec(ctrl = 1, scnr = 1), preset(:decon_2xco2; processes = (humidity = :uniform,));
                        jld2_dir = DATA_DIR, fields = load_climatology(DATA_DIR; dataset = :ncep))
        end
        @test all(r -> all(isfinite, r.Ts), result.ctrl)
        @test isapprox(area_mean(result.scnr), 0.675; atol = 1e-2)
    end
end

@testset "orbital tables: the default rows are near modern; the rcp85 CO2 table loads" begin
    if !isdir(DATA_DIR)
        @test_skip "greb_input_data/ not present"
    else
        modern = quiet(() -> load_climatology(DATA_DIR; dataset = :ncep)).sw_solar
        rms(a) = sqrt(sum(abs2, a .- modern) / length(a))
        for p in (:eccentricity, :obliquity)
            @test rms(resolve(preset(p); jld2_dir = DATA_DIR).solar_table) < 10   # W/m2; row 0 is over 200
        end
        table = resolve(preset(:rcp85); jld2_dir = DATA_DIR).co2_table
        @test isapprox(table[2100], 1231.45; atol = 0.01)
    end
end
