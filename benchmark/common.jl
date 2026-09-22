# Helpers shared by the benchmark scripts.

const REPO = normpath(joinpath(@__DIR__, ".."))
const JULIA_BIN = joinpath(Sys.BINDIR, Base.julia_exename())

"The local dataset directory; never downloads."
function default_data_dir()
    resolved = try
        greb_data_dir(; allow_download=false)
    catch err
        @warn "greb_data_dir could not resolve a dataset; falling back to the repo-local path" err
        nothing
    end
    something(resolved, joinpath(REPO, "greb_input_data"))
end

"Parse a repetition count (an integer >= 1)."
function parse_reps(s::AbstractString)
    n = tryparse(Int, s)
    n === nothing && error("reps must be an integer, got $(repr(s))")
    n >= 1 || error("reps must be at least 1, got $n")
    n
end

"Parse a non-negative integer flag value, e.g. for `--ctrl=10`."
function parse_nonneg_int(name::AbstractString, s::AbstractString)
    n = tryparse(Int, s)
    n === nothing && error("$name must be an integer, got $(repr(s))")
    n >= 0 || error("$name must be >= 0, got $n")
    n
end

"""
Split `--name=value` flags out of `args`. Returns `(flags, positional)`, a
`Dict{String,String}` of flag name to value and the remaining args in order.
Flags can appear anywhere; they don't have to come first.
"""
function split_flags(args::Vector{String})
    flags = Dict{String,String}()
    positional = String[]
    for a in args
        m = match(r"^--([\w-]+)=(.*)$", a)
        if m === nothing
            push!(positional, a)
        else
            flags[m.captures[1]] = m.captures[2]
        end
    end
    return flags, positional
end
