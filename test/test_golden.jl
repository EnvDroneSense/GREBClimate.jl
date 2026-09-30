# Golden regression against a saved snapshot of a real-dataset run.

@testset "golden regression: real dataset control+scenario run matches snapshot" begin
    # Tripwire for any refactor touching the physics kernels: a real 1yr
    # control + 1yr scenario run against the actual NCEP dataset,
    # snapshotted as monthly global-mean Ts/Ta/q, plus a flux = 1 control
    # year so the flux-correction spin-up is exercised too. The tolerances
    # (1e-3 K, 1e-6 kg/kg) are about 30 times the drift measured between runs;
    # they catch behaviour changes, not bit-level drift. Exact equality is
    # checked by tools/validation/bit_identity.jl.
    # Set RUN_GOLDEN=0 to skip this locally. CI skips it too - it has no
    # dataset, so the !isdir(DATA_DIR) branch below always fires there. This
    # guards nothing in CI: a golden break is local-red and CI-green.
    if !isdir(DATA_DIR)
        @test_skip "greb_input_data/ not present"
    elseif get(ENV, "RUN_GOLDEN", "1") == "0"
        @test_skip "RUN_GOLDEN=0"
    else
        fields = load_greb_jld2!(DATA_DIR; dataset = :ncep)
        cfg = create_experiment_config(:full_model)
        result = quiet() do
            greb_model!(RunSpec(flux = 0), cfg; jld2_dir = DATA_DIR, fields = fields)
        end

        gmean(x) = sum(x) / length(x)
        summarize(rec) = (Ts = gmean(rec.Ts), Ta = gmean(rec.Ta), q = gmean(rec.q))

        ctrl_ref = [
            (Ts = 276.6376, Ta = 279.01324, q = 0.006483287),
            (Ts = 276.08505, Ta = 278.59598, q = 0.0066876207),
            (Ts = 275.80753, Ta = 278.15857, q = 0.0068251393),
            (Ts = 276.7514, Ta = 278.90686, q = 0.0069842925),
            (Ts = 278.5325, Ta = 280.68484, q = 0.007306205),
            (Ts = 280.12466, Ta = 282.42026, q = 0.007825737),
            (Ts = 280.6912, Ta = 283.14557, q = 0.008249164),
            (Ts = 280.36264, Ta = 282.90518, q = 0.008242705),
            (Ts = 279.22748, Ta = 281.76135, q = 0.007820126),
            (Ts = 278.26984, Ta = 280.7673, q = 0.0074224365),
            (Ts = 277.96216, Ta = 280.53217, q = 0.0072685555),
            (Ts = 277.9465, Ta = 280.61523, q = 0.007346397),
        ]
        scnr_ref = [
            (Ts = -0.0076227454, Ta = -0.0071252817, q = -4.826931e-08),
            (Ts = -0.0013025337, Ta = -0.0015087194, q = 1.1220921e-08),
            (Ts = -0.00011379851, Ta = -0.0001452234, q = -7.757092e-09),
            (Ts = 0.000121321944, Ta = 0.000109407636, q = 1.9124322e-09),
            (Ts = 0.00016302532, Ta = 0.00015985966, q = 1.4248567e-08),
            (Ts = 9.6678734e-05, Ta = 0.000101053054, q = 1.7348425e-08),
            (Ts = 5.298853e-05, Ta = 5.4988595e-05, q = 1.2192554e-08),
            (Ts = 3.7090645e-05, Ta = 3.6418438e-05, q = 5.587703e-09),
            (Ts = 8.727444e-05, Ta = 8.031726e-05, q = 5.9890044e-09),
            (Ts = 0.00010110272, Ta = 0.00010247363, q = 2.6633937e-09),
            (Ts = 7.9893405e-05, Ta = 8.220805e-05, q = -1.2187116e-09),
            (Ts = 5.7422454e-05, Ta = 5.850858e-05, q = -2.5535258e-09),
        ]

        @test length(result.ctrl) == length(ctrl_ref)
        @test length(result.scnr) == length(scnr_ref)
        for (rec, ref) in zip(result.ctrl, ctrl_ref)
            s = summarize(rec)
            @test isapprox(s.Ts, ref.Ts; atol = 1e-3)
            @test isapprox(s.Ta, ref.Ta; atol = 1e-3)
            @test isapprox(s.q, ref.q; atol = 1e-6)
        end
        for (rec, ref) in zip(result.scnr, scnr_ref)
            s = summarize(rec)
            @test isapprox(s.Ts, ref.Ts; atol = 1e-3)
            @test isapprox(s.Ta, ref.Ta; atol = 1e-3)
            @test isapprox(s.q, ref.q; atol = 1e-6)
        end

        # One spin-up year, then the control: exercises qflux_correction!.
        # Reuses `fields`, which greb_model! restores after the run above.
        flux_result = quiet() do
            greb_model!(RunSpec(flux = 1, ctrl = 1, scnr = 0), cfg; jld2_dir = DATA_DIR, fields = fields)
        end
        flux_ref = [
            (Ts = 277.28638, Ta = 279.7823, q = 0.0063898247),
            (Ts = 276.35574, Ta = 278.86285, q = 0.0064106397),
            (Ts = 275.92465, Ta = 278.2391, q = 0.006436156),
            (Ts = 276.80832, Ta = 278.9109, q = 0.0065385858),
            (Ts = 278.5104, Ta = 280.59595, q = 0.006805867),
            (Ts = 279.90192, Ta = 282.09723, q = 0.0072038374),
            (Ts = 280.28174, Ta = 282.58386, q = 0.007454401),
            (Ts = 279.82877, Ta = 282.20816, q = 0.007350458),
            (Ts = 278.68015, Ta = 281.04367, q = 0.0068896553),
            (Ts = 277.73822, Ta = 280.0452, q = 0.0064565144),
            (Ts = 277.4416, Ta = 279.80316, q = 0.006277522),
            (Ts = 277.36218, Ta = 279.8118, q = 0.0062873242),
        ]
        @test length(flux_result.ctrl) == length(flux_ref)
        for (rec, ref) in zip(flux_result.ctrl, flux_ref)
            s = summarize(rec)
            @test isapprox(s.Ts, ref.Ts; atol = 1e-3)
            @test isapprox(s.Ta, ref.Ta; atol = 1e-3)
            @test isapprox(s.q, ref.q; atol = 1e-6)
        end
    end
end
