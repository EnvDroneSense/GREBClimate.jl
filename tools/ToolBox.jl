# Shared parts of the maintainer scripts: the local dataset, a quiet run, and
# observers that keep one result per phase or sum area-weighted flow terms.
#
# A script loads it with
#   include(joinpath(@__DIR__, "..", "ToolBox.jl")); using .ToolBox
#
# An observer is the `observer(point, view)` function of `greb_model!`; the
# core takes only one, so `both` calls two.

module ToolBox

using GREBClimate
using Logging: with_logger, NullLogger

export data_dir, quiet, repo_state, parse_args, intflag, floatflag, cell_mean, per_phase, both, TermMeans, means

const G = GREBClimate
const REPO = normpath(joinpath(@__DIR__, ".."))

"The local dataset folder; an error when there is none. Never downloads."
function data_dir()
    dir = greb_data_dir(; allow_download=false)
    dir === nothing && error("no local dataset found")
    return dir
end

"Run `f()` with log messages thrown away."
quiet(f) = with_logger(f, NullLogger())

"Short commit hash of the repository and whether tracked files are modified."
function repo_state()
    try
        commit = readchomp(`git -C $REPO rev-parse --short HEAD`)
        dirty = !isempty(readchomp(`git -C $REPO status --porcelain --untracked-files=no`))
        return commit, dirty
    catch
        return "unknown", false
    end
end

"""
    parse_args(args, flags, default) -> (mode, opts, positional)

Split a command line `[mode] [--name=value ...] [positional ...]`. `flags` maps each
mode to the flag names it takes; a first argument that is not a mode, or is a flag,
selects `default`. A bare `--name` gets the value `"true"`. A flag the mode does not
take, or an unknown mode, is an error.
"""
function parse_args(args::Vector{String}, flags::AbstractDict{String}, default::String)
    has_mode = !(isempty(args) || startswith(args[1], "--"))
    mode = has_mode ? args[1] : default
    haskey(flags, mode) || error("unknown mode $(repr(mode)); expected one of $(join(sort!(collect(keys(flags))), ", "))")
    opts = Dict{String,String}()
    positional = String[]
    for a in (has_mode ? args[2:end] : args)
        m = match(r"^--([\w-]+)(?:=(.*))?$", a)
        if m === nothing
            push!(positional, a)
        else
            m.captures[1] in flags[mode] || error("the $mode mode takes no --$(m.captures[1]) flag")
            opts[m.captures[1]] = something(m.captures[2], "true")
        end
    end
    return mode, opts, positional
end

"The integer flag `name` of at least `min`, or `default` when it is absent."
function intflag(opts::AbstractDict, name::AbstractString, default; min::Int=0)
    haskey(opts, name) || return default
    n = tryparse(Int, opts[name])
    n === nothing && error("--$name must be an integer, got $(repr(opts[name]))")
    n >= min || error("--$name must be at least $min, got $n")
    return n
end

"The positive number flag `name`, or `default` when it is absent."
function floatflag(opts::AbstractDict, name::AbstractString, default)
    haskey(opts, name) || return default
    x = tryparse(Float64, opts[name])
    (x === nothing || x <= 0) && error("--$name must be a positive number, got $(repr(opts[name]))")
    return x
end

const WEIGHT = [G.dxlat_grid[j] for _ in 1:G.xdim, j in 1:G.ydim] ./ (G.xdim * sum(G.dxlat_grid))

"Area-weighted global mean of `f(i, j)` over the grid. A per-cell closure keeps the terms free of temporary arrays."
cell_mean(f) = sum(f(i, j) * WEIGHT[i, j] for j in 1:G.ydim, i in 1:G.xdim)

"Observer that creates one `make()` per phase on first use and forwards each call to it. Read one back with `p[phase]`."
struct PerPhase{F,T}
    make::F
    by_phase::Dict{Any,T}
end
per_phase(make) = PerPhase(make, Dict{Any,typeof(make())}())
(p::PerPhase)(point, view) = get!(p.make, p.by_phase, view.phase)(point, view)
Base.getindex(p::PerPhase, phase) = p.by_phase[phase]

"Observer that calls `a`, then `b`. The two stay in the fields `a` and `b`."
struct Both{A,B}
    a::A
    b::B
end
both(a, b) = Both(a, b)
function (o::Both)(point, view)
    o.a(point, view)
    o.b(point, view)
    return nothing
end

"""
Observer that sums `term(view, cap_surf)` of every `terms` entry (`label => term`)
at `at`. `means(t)` divides the sums by the number of steps seen. `cap_surf` is
read before `seaice!` changes it, as the update uses it.
"""
mutable struct TermMeans
    terms::Vector{Pair{String,Function}}
    at::Symbol
    sums::Vector{Float64}
    steps::Int
end
TermMeans(terms; at=:after_tendencies) = TermMeans(terms, at, zeros(length(terms)), 0)
function (t::TermMeans)(point, view)
    point === t.at || return nothing
    for (k, (_, term)) in enumerate(t.terms)
        t.sums[k] += term(view, view.fields.cap_surf)
    end
    t.steps += 1
    return nothing
end
means(t::TermMeans) = t.sums ./ max(t.steps, 1)

end # module
