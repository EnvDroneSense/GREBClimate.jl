### Global-mean shortwave change per unit stratospheric optical depth ###
#
# MAINTAINER TOOL - not part of the package. Needs the local dataset.
#
# Applies a latitude-uniform optical depth through the delta-Eddington
# transmission (the defaults of `AerosolScenario`) and reports the change in
# annual, global-mean net shortwave from `SWradiation!`, evaluated on the
# surface-temperature climatology at every timestep of the year, once
# uncalibrated (`tau_scale` 1) and once with the default `tau_scale`. The result
# is compared with the literature value of about -30 W/m2 per unit optical depth
# (Sato et al. 1993, citing Lacis et al. 1992); see the calibration section of
# the aerosol kernel spec in the vault.
#
# Usage:
#   julia --project=. tools/diagnostics/aerosol_forcing_per_aod.jl

using GREBClimate
using Printf

const REFERENCE = -30.0      # W/m2 per unit optical depth

include(joinpath(@__DIR__, "..", "common.jl"))

# Annual, area-weighted global means of the net shortwave, the incoming
# shortwave and the combined albedo, with the multiplier set to `mult`.
function global_means(fields, cfg, mult)
    state = ModelState()
    state.sw_solar_forcing .= mult
    ws = CirculationWorkspace()
    ts = TimeState(1, 1)
    w = cosd.(Float64.(GREBClimate.lat_grid))
    sw = incoming = albedo = 0.0
    for k in 1:nstep_yr
        ts.ityr = k
        out = SWradiation!(@view(fields.Tclim[:, :, k]), fields, state, ts, cfg, ws)
        for j in 1:ydim
            s0 = Float64(fields.sw_solar[j, k]) * mult[j] * 0.01 * GREBClimate.S0_var
            for i in 1:xdim
                sw += w[j] * out.SW[i, j]
                incoming += w[j] * s0
                albedo += w[j] * out.albedo[i, j]
            end
        end
    end
    norm = nstep_yr * xdim * sum(w)
    return (sw=sw / norm, incoming=incoming / norm, albedo=albedo / norm)
end

function main()
    fields = load_greb_jld2!(DATA_DIR; dataset=:ncep)
    cfg = create_experiment_config(:full_model)
    sc = AerosolScenario()
    base = global_means(fields, cfg, ones(Float32, ydim))
    @printf("control: net SW %.1f W/m2, incoming %.1f W/m2, albedo %.3f\n",
            base.sw, base.incoming, base.albedo)
    @printf("transmission: ssa %.2f, g %.2f, mu0 %.2f\n", sc.ssa, sc.asymmetry, sc.mu0)
    for scale in (1.0, sc.tau_scale)
        @printf("\ntau_scale %.3f%s\n", scale, scale == 1.0 ? " (uncalibrated)" : " (default)")
        println("optical depth  R/tau   dSW [W/m2]  dSW/tau [W/m2]  ratio to reference")
        for tau in (0.01, 0.05, 0.1, 0.2, 0.5)
            R, T = GREBClimate.delta_eddington(scale * tau, sc.ssa, sc.asymmetry, sc.mu0)
            d = global_means(fields, cfg, fill(Float32(T), ydim)).sw - base.sw
            @printf("%-13.2f  %.3f   %-10.3f  %-14.1f  %.2f\n",
                    tau, R / tau, d, d / tau, d / tau / REFERENCE)
        end
    end
end

main()
