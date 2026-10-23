# The names the tests call besides the exports. Shared by runtests.jl, live.jl
# and the subprocesses of test_threading.jl.
using GREBClimate
using GREBClimate:
    # The kernels and the run loop are internal, reached by name
    SWradiation!, LWradiation!, hydro!, convergence!, seaice!, deep_ocean!,
    diffusion!, advection!, circulation!, tendencies!, time_loop!, output!, diagnostics!,
    qflux_correction!, init_model!, derive_fields!, _mutated_fields,
    # The calendar
    spinup, control, scenario, steps, year, month, months_per_year, month_of_day,
    is_day_end, is_month_end, is_year_end, steps_in_month, first_step_of,
    data_slice, previous_slice, clim_slice, RecordTime,
    # State, output and checks
    OUTPUT_FIELDS, BoundaryAnomalyFields, BudgetCheck, RangeCheck, in_range, is_land, jldopen,
    # The grid and the physical constants the kernel tests calculate with
    lat_grid, dlon, z_air, z_vapor, min_T_K, cap_land, cap_ocean, d_ocean,
    To_ice1, To_ice2, Tl_ice1, Tl_ice2, inv_To_ice_range, inv_Tl_ice_range,
    a_no_ice, da_ice, a_cloud, emissivity_fit, σ, ρ_air, ce, cq_latent, cq_rain, r_qviwv,
    ct_sens, turb_coeff, c_effmix, q_to_mm_per_day, q_to_mm_per_day
using Logging: with_logger, NullLogger
using Test
