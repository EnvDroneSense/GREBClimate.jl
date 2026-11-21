# Same-process timing of two small parts of a step:
#   forcing - the per-step CO2 and solar lookup, against an empty call
#   seaice! - the loop order (i outer, as in src/) against a j-outer copy
#
#   julia --project=benchmark -t 2,0 benchmark/compare_step_parts.jl [jld2_dir] [--seconds=S]
#
# Run it on a quiet machine, ideally pinned (benchmark/README.md). It stops if
# the two seaice! variants differ. A ratio within about 2 % of 1 is noise.

using GREBClimate
using GREBClimate: seaice!, forcing, ModelTime, scenario, xdim, ydim, clim_slice, cap_ocean, cap_land,
                   To_ice1, To_ice2, inv_To_ice_range, @is_land, @ice_ramp

include("common.jl")
include("compare.jl")

# `seaice!` with the ice-albedo feedback, looping j outer and i inner (column-major order)
function seaice_j_outer!(fields, Ts, slice)
    mld = clim_slice(fields.mld_clim, slice)
    cap_surf = fields.cap_surf
    @inbounds @fastmath for j in 1:ydim, i in 1:xdim

        ice_frac = @ice_ramp(Ts[i, j], To_ice1, To_ice2, inv_To_ice_range)
        cap_ice = cap_land * ice_frac + cap_ocean * mld[i, j] * (1.0f0 - ice_frac)
        cap_surf[i, j] = ifelse(@is_land(fields.z_topo[i, j]), cap_surf[i, j], cap_ice)
    end
    @. cap_surf = ifelse(fields.glacier > 0.5f0, cap_land, cap_surf)
    return nothing
end

function main(args)
    flags, pos = split_flags(args)
    check_flags("compare_step_parts", flags, ("seconds",))
    seconds = parse(Float64, get(flags, "seconds", "5"))
    dir = check_data_dir(isempty(pos) ? default_data_dir() : pos[1])
    (; fields, r, slice, Ts) = stage_inputs(dir, preset(:full_model))
    p = r.config.processes
    t = ModelTime(scenario, 1950, 1)

    seaice!(fields, Ts, slice, p)
    expected = copy(fields.cap_surf)
    seaice_j_outer!(fields, Ts, slice)
    fields.cap_surf == expected || error("seaice_j_outer! differs from seaice!")

    calibration = machine_header()
    ok_seaice = compare_kernels("seaice! (i outer)" => () -> seaice!(fields, Ts, slice, p),
        "seaice! (j outer)" => () -> seaice_j_outer!(fields, Ts, slice); seconds)
    ok_forcing = compare_kernels("forcing" => () -> forcing(t, r), "empty call" => () -> nothing; seconds)
    machine_footer(calibration)
    return ok_seaice && ok_forcing ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
