module GREBClimate

# =============================================================================
# GREB - Globally Resolved Energy Balance model
#
# A global climate model that steps surface/air/ocean temperature and
# humidity forward under shortwave/longwave radiation, hydrology, sea ice,
# deep-ocean coupling, and atmospheric circulation (diffusion, advection,
# moisture convergence). An interactive Pluto version of the same model
# lives in `notebooks/GREB_explorer.jl`; see the package docs or README for
# usage.
#
# Files below are included in dependency order: constants -> calendar ->
# config/{processes,scenario,presets} -> data -> state -> io -> config/resolve
# -> physics/{radiation,hydrology,ocean,circulation} -> tendencies -> budgets
# -> forcing -> output -> postprocess -> model -> ensemble.
# =============================================================================

using LoopVectorization   # @turbo SIMD
using JLD2
using Logging: with_logger, NullLogger
using DataDeps: DataDeps, DataDep, register, unpack, @datadep_str

export Config, preset, preset_names, Scenario, Processes, Hydrology, mscm_hydrology
export resolve, ResolvedConfig, ResolvedHydrology
export Corrections, SpinUp, Stored, NoCorrections
export CO2Path, ConstantCO2, CO2Table, CO2File, A1BRamp, CO2SineWave, CO2Step, SeasonalCO2
export CO2Mask, UniformMask, LatitudeMask, SurfaceMask
export Solar, ModernSolar, SolarConstant, SolarCycle, SolarTable, EarthSunDistance
export SurfaceForcing, NoSurfaceForcing, BoundaryAnomaly, SSTOffset
export RunSpec, CirculationWorkspace, MonthlyAccumulator, TimeState, MonthlyRecord
export ClimateFields, ModelState, SurfaceState
export greb_data_dir
export load_climatology, load_flux_corrections!, load_boundary_anomaly!, read_field
export load_co2_scenario, load_co2_custom, load_solar_forcing
export apply_co2_mask!, apply_dynamic_co2_mask!, forcing
export build_monthly_climatology, apply_scenario_anomalies, compute_annual_ice_climatology, global_mean
export greb_model!, run_ensemble
export xdim, ydim, nstep_yr
export step_of_year, day_of_year, month_of_step, decimal_year

include("core/constants.jl")
include("core/calendar.jl")
include("config/processes.jl")
include("config/scenario.jl")
include("config/presets.jl")
include("data.jl")
include("core/state.jl")
include("io.jl")
include("config/resolve.jl")
include("physics/radiation.jl")
include("physics/hydrology.jl")
include("physics/ocean.jl")
include("physics/circulation.jl")
include("core/tendencies.jl")
include("core/budgets.jl")
include("forcing/forcing.jl")
include("core/output.jl")
include("core/postprocess.jl")
include("core/model.jl")
include("core/ensemble.jl")

function __init__()
    # Registration only: nothing is downloaded until `greb_data_dir()` has to
    # fall through to the DataDep, which cannot happen during precompilation.
    register_greb_datadep()
    return nothing
end

using PrecompileTools: @compile_workload

@compile_workload begin
    with_logger(NullLogger()) do
        # A control and a scenario year, so the scenario phase is compiled too
        greb_model!(RunSpec(ctrl=1, scnr=1), preset(:co2_double; corrections=NoCorrections());
                    jld2_dir="", allow_uninitialized=true)
        greb_model!(RunSpec(scnr=0), preset(:full_model; hydrology=(evaporation=:skin,),
                                             corrections=NoCorrections());
                    jld2_dir="", allow_uninitialized=true)
    end
end

end # module GREBClimate
