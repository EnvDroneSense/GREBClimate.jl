# Golden regression against a saved snapshot of a real-dataset run.

include(joinpath("data", "golden_reference.jl"))

"The monthly global means of `records` match the snapshot `reference`, record by record."
function matches_reference(records, reference)
    @test length(records) == length(reference)
    for (rec, ref) in zip(records, reference)
        @test isapprox(gmean(rec.Ts), ref.Ts; atol=1e-3)
        @test isapprox(gmean(rec.Ta), ref.Ta; atol=1e-3)
        @test isapprox(gmean(rec.q), ref.q; atol=1e-6)
    end
end

@testset "golden regression: real dataset control+scenario run matches snapshot" begin
    # A 1-year control and 1-year scenario on the NCEP dataset, stored as
    # monthly global-mean Ts/Ta/q, plus a control year after a one-year spin-up.
    # The tolerances (1e-3 K, 1e-6 kg/kg) are about 30 times the drift measured
    # between runs; exact equality is checked by tools/validation/bit_identity.jl.
    # RUN_GOLDEN=0 skips it locally. CI has no dataset and always skips it, so
    # a golden break is local-red and CI-green.
    with_dataset() do
        get(ENV, "RUN_GOLDEN", "1") == "0" && return @test_skip "RUN_GOLDEN=0"
        fields = load_climatology(DATA_DIR; dataset=:ncep)
        result = quiet() do
            greb_model!(RunSpec(), preset(:full_model; corrections=Stored()); jld2_dir=DATA_DIR, fields=fields)
        end
        matches_reference(result.ctrl, GOLDEN_CTRL)
        matches_reference(result.scnr, GOLDEN_SCNR)

        # One spin-up year, then the control: exercises qflux_correction!.
        # Reuses `fields`, which greb_model! restores after the run above.
        flux_result = quiet() do
            greb_model!(RunSpec(ctrl=1, scnr=0), preset(:full_model; corrections=SpinUp(1));
                jld2_dir=DATA_DIR, fields=fields)
        end
        matches_reference(flux_result.ctrl, GOLDEN_FLUX)
    end
end

# Area-weighted global mean of a monthly-record field, averaged over `recs`
function area_mean(recs, var=:Ts)
    w = cosd.(Float64.(lat_grid))
    return sum(sum(getfield(r, var) .* w') for r in recs) / (length(recs) * X * sum(w))
end

@testset "MSCM configuration reproduces the MSCM 2xCO2 response" begin
    # MSCM (Monash Simple Climate Model) database, 2xCO2 with every process on:
    # year-1 global-mean surface temperature response 0.594636 K.
    with_dataset() do
        cfg = preset(:co2_double; processes=(moisture_convergence=false,), hydrology=mscm_hydrology())
        result = quiet() do
            greb_model!(RunSpec(ctrl=1, scnr=1), cfg; jld2_dir=DATA_DIR,
                fields=load_climatology(DATA_DIR; dataset=:ncep))
        end
        @test isapprox(area_mean(result.scnr), MSCM_YEAR1; atol=1e-3)
    end
end

@testset "mean-climate deconstruction: a switched-off process changes the control climate" begin
    # On computed corrections every configuration is pulled back to the
    # observed climate; on the stored ones the switch shows.
    with_dataset() do
        fields = load_climatology(DATA_DIR; dataset=:ncep)
        control(processes) = quiet() do
            greb_model!(RunSpec(ctrl=1, scnr=0), preset(:decon_mean_climate; processes);
                jld2_dir=DATA_DIR, fields)
        end.ctrl
        @test abs(area_mean(control((ocean=:none,))) - area_mean(control((;)))) > 0.1
    end
end

@testset "2xCO2 deconstruction: uniform humidity gives a finite response" begin
    # 0.0052 kg/kg is above saturation over cold, high ground. The fitted rain
    # scheme rains that out within the spin-up and the run goes non-finite;
    # the preset's MSCM physics does not.
    with_dataset() do
        result = quiet() do
            greb_model!(RunSpec(ctrl=1, scnr=1), preset(:decon_2xco2; processes=(humidity=:uniform,));
                jld2_dir=DATA_DIR, fields=load_climatology(DATA_DIR; dataset=:ncep))
        end
        @test all(r -> all(isfinite, r.Ts), result.ctrl)
        @test isapprox(area_mean(result.scnr), 0.675; atol=1e-2)
    end
end

@testset "orbital tables: the default rows are near modern; the rcp85 CO2 table loads" begin
    with_dataset() do
        modern = quiet(() -> load_climatology(DATA_DIR; dataset=:ncep)).sw_solar
        rms(a) = sqrt(sum(abs2, a .- modern) / length(a))
        for p in (:eccentricity, :obliquity)
            @test rms(resolve(preset(p); jld2_dir=DATA_DIR).solar_table) < 10   # W/m2; row 0 is over 200
        end
        table = resolve(preset(:rcp85); jld2_dir=DATA_DIR).co2_table
        @test isapprox(table[2100], 1231.45; atol=0.01)
    end
end
