# Energy and water bookkeeping of one run, from the per-step observer.
#
#   julia --project=. -t 2,0 tools/checks/budget.jl [preset] [years]
#
# Runs `preset` (default full_model) on the stored flux corrections for `years`
# (default 1) of control and scenario. Per phase it prints the largest gap
# between each store's change and the sum of its flows, the cell-steps held by a
# limiter (GREBClimate.BudgetCheck), and the global-mean flows. Needs the local
# dataset.

using GREBClimate
using Printf
include(joinpath(@__DIR__, "..", "ToolBox.jl"))
using .ToolBox
const G = GREBClimate

const TERMS = [
    "shortwave absorbed by the surface [W/m2]" => (v, c) -> cell_mean((i, j) -> v.tend.SW[i, j]),
    "longwave out to space (olr) [W/m2]" => (v, c) -> cell_mean((i, j) -> -v.tend.LW_up[i, j] - (1 - v.tend.em[i, j]) * v.tend.LW_surf[i, j]),
    "net flux into the surface [W/m2]" => (v, c) -> cell_mean((i, j) -> G.surface_flux(v.tend, i, j)),
    "surface flux correction [W/m2]" => (v, c) -> cell_mean((i, j) -> v.fields.Ts_flux_correction[i, j, v.step_of_year]),
    "heat from the deep ocean to the surface [W/m2]" => (v, c) -> cell_mean((i, j) -> c[i, j] * v.tend.dT_ocean[i, j] / G.Δt),
    "net flux into the air [W/m2]" => (v, c) -> cell_mean((i, j) -> G.atmosphere_flux(v.tend, i, j)),
    "latent heat: surface loss + air gain [W/m2]" => (v, c) -> cell_mean((i, j) -> v.tend.Q_lat[i, j] + v.tend.Q_lat_air[i, j]),
    "air temperature change by transport [K/step]" => (v, c) -> cell_mean((i, j) -> v.tend.dTa_crcl[i, j]),
    "evaporation minus rain [kg/kg per step]" => (v, c) -> cell_mean((i, j) -> G.Δt * (v.tend.dq_eva[i, j] + v.tend.dq_rain[i, j])),
    "humidity change by transport [kg/kg per step]" => (v, c) -> cell_mean((i, j) -> v.tend.dq_crcl[i, j]),
    "humidity flux correction [kg/kg per step]" => (v, c) -> cell_mean((i, j) -> v.fields.q_flux_correction[i, j, v.step_of_year])
]

function main(args)
    name = Symbol(get(args, 1, "full_model"))
    years = parse(Int, get(args, 2, "1"))
    dir = data_dir()
    fields = load_climatology(dir; dataset=:ncep)
    phases = per_phase(() -> both(TermMeans(TERMS), G.BudgetCheck()))
    quiet() do
        greb_model!(RunSpec(ctrl=years, scnr=years), preset(name; corrections=Stored());
            jld2_dir=dir, fields=fields, observer=phases)
    end
    println("preset :", name, ", ", years, " yr control + ", years, " yr scenario, ",
        Threads.nthreads(), " thread(s)")
    for phase in (G.control, G.scenario)
        terms, b = phases[phase].a, phases[phase].b
        println("\n== ", phase, " (", b.steps, " steps) ==")
        println("largest gap between a store's change and its flows:")
        @printf("  surface %.2e K   air %.2e K   ocean %.2e K   humidity %.2e kg/kg\n",
            b.surface, b.atmosphere, b.ocean, b.water)
        println("limiter cell-steps (of ", b.steps * G.xdim * G.ydim, "):")
        println("  Ts floor ", b.floor_Ts, "   Ta floor ", b.floor_Ta, "   humidity removal ", b.humidity_low,
            "   humidity gain ", b.humidity_high, "   rain limit ", b.rain_limit)
        println("global means:")
        for ((label, _), m) in zip(TERMS, means(terms))
            @printf("  %-48s %+.6g\n", label, m)
        end
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
