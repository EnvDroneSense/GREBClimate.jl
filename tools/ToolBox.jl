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

export data_dir, quiet, cell_mean, per_phase, both, TermMeans, means

const G = GREBClimate

"The local dataset folder; an error when there is none. Never downloads."
function data_dir()
    dir = greb_data_dir(; allow_download=false)
    dir === nothing && error("no local dataset found")
    return dir
end

"Run `f()` with log messages thrown away."
quiet(f) = with_logger(f, NullLogger())

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
