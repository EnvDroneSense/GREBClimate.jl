# State: the CO2 mask reset and the derived fields.

@testset "init_model! clears the CO2 mask an earlier run left in the same fields" begin
    fields = ClimateFields()
    apply_co2_mask!(fields, LatitudeMask(:nh))
    @test any(!=(1.0f0), fields.co2_part)  # what a regional scenario leaves behind

    # Full CO2 again, also for a regional preset: its mask is for the scenario
    for p in (:full_model, :regional_co2_nh)
        fields.co2_part[1, 1] = 0.5f0
        quiet() do
            init_model!(fields, resolve(preset(p)))
        end
        @test all(==(1.0f0), fields.co2_part)
    end
end

@testset "derive_fields! follows the input maps" begin
    fields = synthetic_fields()
    quiet(() -> init_model!(fields, resolve(preset(:full_model))))
    @test fields.cap_surf[60, 10] == cap_ocean * fields.mld_clim[60, 10, 1]

    fields.z_topo[60, 10] = 500.0f0          # an ocean cell becomes land
    fields.mld_clim[70, 20, :] .= 80.0f0
    fields.u_clim[1, 1, 1] = -3.0f0
    fields.v_clim[2, 2, 2] = 4.0f0
    derive_fields!(fields, Processes())

    @test fields.cap_surf[60, 10] == cap_land
    @test fields.wz_air[60, 10] == exp(-500.0f0 / z_air)
    @test fields.wz_vapor[60, 10] == exp(-500.0f0 / z_vapor)
    @test fields.z_ocean[70, 20] == 240.0f0
    @test fields.cap_surf[70, 20] == cap_ocean * 80.0f0
    @test (fields.u_clim_neg[1, 1, 1], fields.u_clim_pos[1, 1, 1]) == (-3.0f0, 0.0f0)
    @test (fields.v_clim_neg[2, 2, 2], fields.v_clim_pos[2, 2, 2]) == (0.0f0, 4.0f0)
    @test fields.u_clim_pos .+ fields.u_clim_neg == fields.u_clim
    @test fields.v_clim_pos .+ fields.v_clim_neg == fields.v_clim
    @test all(>=(0), fields.u_clim_pos) && all(<=(0), fields.u_clim_neg)
end
