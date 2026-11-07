### Bit-identity check for refactors ###
#
# MAINTAINER TOOL. Needs the local dataset.
#
# Saves every monthly record field of the control and scenario runs of every
# experiment preset, plus two cases with switches off, then compares a later
# run with exact equality (NaN equals NaN). A refactor that claims "no change
# to results" must report 0 differing values. A snapshot is about 250 MB; a run
# takes a few minutes.
#
# Snapshots compare only on the same machine, Julia version and thread count
# (`@turbo` differs between CPUs and releases; a thread mismatch is refused).
# Save on the commit before the change, compare on the commit after, once with
# `-t 1` and once with `-t 2,0`.
#
# Usage:
#   julia --project=. tools/validation/bit_identity.jl save    <snapshot.jld2>
#   julia --project=. tools/validation/bit_identity.jl compare <snapshot.jld2>

using GREBClimate
using Logging: with_logger, NullLogger

function data_dir()
    dir = greb_data_dir(; allow_download=false)
    dir === nothing && error("no local dataset found")
    return dir
end

# One spin-up year (exercises the flux-correction loop), one control and one
# scenario year. A snapshot must hold every case of the run: retake it when a
# preset is renamed.
const RUN = RunSpec(ctrl=1, scnr=1)

case(name; kw...) = preset(name; corrections=SpinUp(1), kw...)

# CO2 path for the :custom_co2 preset
function _custom_co2_file()
    path = joinpath(mktempdir(), "custom_co2.txt")
    write(path, "1950 400\n1951 420\n")
    return path
end

# name => config constructor: every preset, then the two switch cases
function _cases()
    cases = Pair{String,Function}[]
    for p in preset_names()
        mk = p === :custom_co2 ? () -> case(p; path=_custom_co2_file()) : () -> case(p)
        push!(cases, string(p) => mk)
    end
    push!(cases, "flat_topography" => () -> preset(:co2_double; processes=(topography=:flat,), corrections=Stored()))
    push!(cases, "decon_crcl_hydro_off" => () -> case(:full_model; processes=(transport=false, hydrology=:none)))
    return sort!(cases; by=first)
end

function run_cases()
    dir = data_dir()
    fields = load_climatology(dir; dataset=:ncep)
    snap = Dict{String,Array{Float32,3}}()
    for (name, mkcfg) in _cases()
        result = with_logger(NullLogger()) do
            greb_model!(RUN, mkcfg(); jld2_dir=dir, fields=fields)
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
    jv == string(VERSION) || @warn "snapshot from Julia $jv; exact equality is not expected across versions"
    new = run_cases()
    missing_now = setdiff(keys(old), keys(new))
    isempty(missing_now) || error("the run lacks $(length(missing_now)) snapshot arrays, e.g. $(first(missing_now))")
    case_of(k) = first(split(k, '/'))
    old_cases = Set(case_of(k) for k in keys(old))
    added = setdiff(keys(new), keys(old))
    new_cases = sort!(unique(c for c in (case_of(k) for k in added) if !(c in old_cases)))
    new_fields = sort!(unique(last(split(k, '/')) for k in added if case_of(k) in old_cases))
    isempty(new_cases) || println("  new cases (not compared): ", join(new_cases, ", "))
    isempty(new_fields) || println("  new record fields (not compared): ", join(new_fields, ", "))
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

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
