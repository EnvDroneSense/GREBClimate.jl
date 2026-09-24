### GREB's temperature response to a uniform shortwave forcing ###
#
# MAINTAINER TOOL - not part of the package. Needs the local dataset.
#
# Dims the shortwave by a uniform factor through `cfg.solar_scenario`, chosen so
# the global-mean absorbed shortwave falls by FORCING W/m2, once for one year
# (a volcano-like pulse) and once from that year on (a step). Reports the
# global-mean surface temperature response per W/m2, its timing and recovery,
# land and ocean separately, and the mean mixed-layer depth that sets the
# ocean's heat capacity. Separates GREB's response to a known forcing from the
# forcing the aerosol add-on produces; see the vault note volcanic-overcooling.
#
# Usage:
#   julia --project=. tools/shortwave_pulse_response.jl

using GREBClimate
using Printf

const REPO = normpath(joinpath(@__DIR__, ".."))
const DATA_DIR = something(greb_data_dir(; allow_download=false),
                           joinpath(REPO, "greb_input_data"))
const FORCING = 3.0          # W/m2, about the Pinatubo peak effective forcing
const NET_SW = 233.8         # control global-mean absorbed shortwave (tools/aerosol_forcing_per_aod.jl)
const YEARS = 30
const START = 1951           # forcing year; the scenario starts in 1950

w_lat() = cosd.(Float64.(GREBClimate.lat_grid))

function means(result, mask)
    w = reshape(w_lat(), 1, :) .* mask
    return [sum(rec.Ts .* w) / sum(w) for rec in result.scnr]
end

function run(fields, factor_of_year)
    cfg = create_experiment_config(:full_model)
    cfg.solar_scenario = Dict(y => Float32(factor_of_year(y)) for y in 1950:1950 + YEARS)
    return redirect_stdout(devnull) do
        greb_model!(RunSpec(flux=3, ctrl=5, scnr=YEARS), cfg; jld2_dir=DATA_DIR,
                    fields=deepcopy(fields))
    end
end

function main()
    fields = load_greb_jld2!(DATA_DIR; dataset=:ncep)
    land = Float64.(fields.z_topo .> 0)
    masks = (("global", ones(xdim, ydim)), ("land", land), ("ocean", 1 .- land))
    x = 1 - FORCING / NET_SW
    base = run(fields, y -> 1.0)
    pulse = run(fields, y -> y == START ? x : 1.0)
    step = run(fields, y -> y >= START ? x : 1.0)
    m0 = 12 * (START - 1950)                      # last month before the forcing

    @printf("uniform shortwave factor %.5f (-%.1f W/m2 global-mean absorbed)\n", x, FORCING)
    for (name, mask) in masks
        b = means(base, mask)
        p = means(pulse, mask) .- b
        s = means(step, mask) .- b
        i = argmin(p)
        rec = findfirst(v -> v > p[i] / exp(1), p[i:end])
        @printf("\n%s\n", name)
        @printf("  one-year pulse: peak %.3f K (%.3f K per W/m2) in month %d; 1/e recovery %s\n",
                p[i], p[i] / FORCING, i - m0,
                rec === nothing ? "not within the run" : "$(rec - 1) months after the peak")
        @printf("  pulse tail, annual mean: %s K after 2, 5 and 10 years\n",
                join((@sprintf("%.3f", sum(p[m0+12k+1:m0+12k+12]) / 12) for k in (2, 5, 10)), ", "))
        for yrs in (1, 2, 5, 10, 20)
            m = m0 + 12 * yrs
            m <= length(s) || continue
            @printf("  step, after %2d yr: %.3f K (%.3f K per W/m2)\n", yrs,
                    sum(s[m-11:m]) / 12, sum(s[m-11:m]) / 12 / FORCING)
        end
    end
    ocean = fields.z_topo .<= 0
    w = reshape(w_lat(), 1, :) .* ocean
    mld = sum(sum(fields.mldclim; dims=3)[:, :, 1] ./ nstep_yr .* w) / sum(w)
    @printf("\narea-mean annual mixed-layer depth over ocean: %.1f m\n", mld)
end

main()
