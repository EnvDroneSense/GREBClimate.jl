# Energy and water bookkeeping of one run, from the per-step observer.
#
#   julia --project=. -t 2,0 tools/diagnostics/budget.jl [preset] [years]
#
# Runs `preset` (default full_model) on the stored flux corrections for
# `years` (default 1) of control and scenario and prints, per phase:
#   - the largest difference between each store's change and the sum of its
#     flows, and how many cells a limiter held (GREBClimate.BudgetCheck);
#   - the global-mean flows, averaged over the phase.
# Needs the local dataset.

using GREBClimate
using Printf
const G = GREBClimate

const WEIGHT = [G.dxlat_grid[j] for _ in 1:G.xdim, j in 1:G.ydim] ./ (G.xdim * sum(G.dxlat_grid))
gmean(f) = sum(f(i, j) * WEIGHT[i, j] for j in 1:G.ydim, i in 1:G.xdim)

const TERMS = [
    "shortwave absorbed by the surface [W/m2]"      => (v, c) -> gmean((i, j) -> v.tend.SW[i, j]),
    "longwave out to space (olr) [W/m2]"            => (v, c) -> gmean((i, j) -> -v.tend.LW_up[i, j] - (1 - v.tend.em[i, j]) * v.tend.LW_surf[i, j]),
    "net flux into the surface [W/m2]"              => (v, c) -> gmean((i, j) -> G.surface_flux(v.tend, i, j)),
    "surface flux correction [W/m2]"                => (v, c) -> gmean((i, j) -> v.fields.Ts_flux_correction[i, j, v.step_of_year]),
    "heat from the deep ocean to the surface [W/m2]" => (v, c) -> gmean((i, j) -> c[i, j] * v.tend.dT_ocean[i, j] / G.Δt),
    "net flux into the air [W/m2]"                  => (v, c) -> gmean((i, j) -> G.atmosphere_flux(v.tend, i, j)),
    "latent heat: surface loss + air gain [W/m2]"   => (v, c) -> gmean((i, j) -> v.tend.Q_lat[i, j] + v.tend.Q_lat_air[i, j]),
    "air temperature change by transport [K/step]"  => (v, c) -> gmean((i, j) -> v.tend.dTa_crcl[i, j]),
    "evaporation minus rain [kg/kg per step]"       => (v, c) -> gmean((i, j) -> G.Δt * (v.tend.dq_eva[i, j] + v.tend.dq_rain[i, j])),
    "humidity change by transport [kg/kg per step]" => (v, c) -> gmean((i, j) -> v.tend.dq_crcl[i, j]),
    "humidity flux correction [kg/kg per step]"     => (v, c) -> gmean((i, j) -> v.fields.q_flux_correction[i, j, v.step_of_year]),
]

mutable struct Phase
    check::G.BudgetCheck
    sums::Vector{Float64}
end
Phase() = Phase(G.BudgetCheck(), zeros(length(TERMS)))

function main(args)
    name = Symbol(get(args, 1, "full_model"))
    years = parse(Int, get(args, 2, "1"))
    dir = greb_data_dir(; allow_download=false)
    dir === nothing && error("no local dataset found")
    fields = load_climatology(dir; dataset=:ncep)
    phases = Dict(G.control => Phase(), G.scenario => Phase())
    function observer(point, view)
        ph = phases[view.phase]
        if point === :after_tendencies
            # cap_surf as the update uses it; seaice! changes it afterwards
            for (k, (_, term)) in enumerate(TERMS)
                ph.sums[k] += term(view, view.fields.cap_surf)
            end
        end
        ph.check(point, view)
    end
    Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        greb_model!(RunSpec(ctrl=years, scnr=years), preset(name; corrections=Stored());
            jld2_dir=dir, fields=fields, observer=observer)
    end
    println("preset :", name, ", ", years, " yr control + ", years, " yr scenario, ",
        Threads.nthreads(), " thread(s)")
    for phase in (G.control, G.scenario)
        ph = phases[phase]; b = ph.check
        println("\n== ", phase, " (", b.steps, " steps) ==")
        println("largest difference between a store's change and its flows, any cell and step:")
        @printf("  surface %.2e K   air %.2e K   ocean %.2e K   humidity %.2e kg/kg\n",
            b.surface, b.atmosphere, b.ocean, b.water)
        println("cells held by a limiter (cell-steps of ", b.steps * G.xdim * G.ydim, "):")
        println("  Ts floor ", b.floor_Ts, "   Ta floor ", b.floor_Ta, "   humidity removal ", b.humidity_low,
            "   humidity gain ", b.humidity_high, "   rain limit ", b.rain_limit)
        println("global means over the phase:")
        for (k, (label, _)) in enumerate(TERMS)
            @printf("  %-48s %+.6g\n", label, ph.sums[k] / max(b.steps, 1))
        end
    end
    return 0
end

exit(main(ARGS))
