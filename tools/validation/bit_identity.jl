### Bit-identity check for refactors ###
#
# MAINTAINER TOOL - not part of the package. Needs the local dataset.
#
# Runs every experiment preset, plus one case with switches off, and saves every
# monthly record field of the control and scenario runs; after a change, runs
# them again and compares with exact equality (NaN equals NaN). A refactor that
# claims "no change to results" must report 0 differing values here. A snapshot
# is about 250 MB and a run takes a few minutes.
#
# Snapshots are only comparable on the same machine and Julia version:
# `@turbo` code differs between CPUs (AVX2, AVX-512) and Julia releases. Take
# the snapshot on the commit before the change, compare on the commit after.
# Run once with `-t 1` and once with `-t 2`; the thread count is stored and a
# mismatch is refused.
#
# Usage:
#   julia --project=. tools/validation/bit_identity.jl save    <snapshot.jld2>
#   julia --project=. tools/validation/bit_identity.jl compare <snapshot.jld2>

using GREBClimate

const DATA_DIR = something(greb_data_dir(; allow_download=false),
                           joinpath(@__DIR__, "..", "..", "greb_input_data"))
isdir(DATA_DIR) || error("dataset not found at $DATA_DIR")

# flux = 1 so the flux-correction spin-up loop is exercised, not only the run loop.
const RUN = RunSpec(flux=1, ctrl=1, scnr=1)

function _decon_cfg()
    cfg = create_experiment_config(:full_model)
    cfg.log_crcl_dmc = false   # the humidity-circulation selection in time_loop!
    cfg.log_hydro_dmc = false
    return cfg
end

# A fixed two-year CO2 path for the :custom_co2 preset
function _custom_co2_file()
    path = joinpath(mktempdir(), "custom_co2.txt")
    write(path, "1950 400\n1951 420\n")
    return path
end

# name => config constructor: every preset, then the switch case
function _cases()
    presets = sort!(collect(keys(GREBClimate._EXPERIMENT_OVERRIDES)))
    append!(presets, (:decon_mean_climate, :decon_2xco2))
    cases = Pair{String,Function}[]
    for p in presets
        mk = p === :custom_co2 ? () -> create_experiment_config(p; co2_path=_custom_co2_file()) :
                                 () -> create_experiment_config(p)
        push!(cases, string(p) => mk)
    end
    push!(cases, "decon_crcl_hydro_off" => _decon_cfg)
    return cases
end

function run_cases()
    fields = load_greb_jld2!(DATA_DIR; dataset=:ncep)
    snap = Dict{String,Array{Float32,3}}()
    for (name, mkcfg) in _cases()
        result = redirect_stdout(devnull) do
            greb_model!(RUN, mkcfg(); jld2_dir=DATA_DIR, fields=deepcopy(fields))
        end
        for phase in (:ctrl, :scnr), var in fieldnames(MonthlyRecord)
            recs = getfield(result, phase)
            snap["$name/$phase/$var"] = cat((getfield(r, var) for r in recs)...; dims=3)
        end
        println("  ran $name")
    end
    return snap
end

function main(args)
    length(args) == 2 && args[1] in ("save", "compare") ||
        error("usage: bit_identity.jl save|compare <snapshot.jld2>")
    mode, path = args
    println("Threads.nthreads() = ", Threads.nthreads(), ", Julia ", VERSION)
    if mode == "save"
        snap = run_cases()
        GREBClimate.jldopen(path, "w") do f
            f["nthreads"] = Threads.nthreads()
            f["julia"] = string(VERSION)
            f["data"] = snap
        end
        println("saved $(length(snap)) arrays to $path")
        return 0
    end
    old, nt, jv = GREBClimate.jldopen(path, "r") do f
        f["data"], f["nthreads"], f["julia"]
    end
    nt == Threads.nthreads() || error("snapshot taken with $nt threads; rerun with -t $nt")
    jv == string(VERSION) || @warn "snapshot taken on Julia $jv; exact equality is not expected across versions"
    new = run_cases()
    keys(old) == keys(new) || error("snapshot and current run hold different arrays")
    total = 0
    for k in sort!(collect(keys(old)))
        a, b = old[k], new[k]
        size(a) == size(b) || (println("  $k: size $(size(a)) -> $(size(b))"); total += 1; continue)
        n = count(!isequal(x, y) for (x, y) in zip(a, b))
        if n > 0
            println("  $k: $n differing values, max |diff| $(maximum(abs, a .- b))")
            total += n
        end
    end
    println(total == 0 ? "IDENTICAL: 0 differing values in $(length(old)) arrays" :
                         "DIFFERENT: $total differing values")
    return total == 0 ? 0 : 1
end

exit(main(ARGS))
