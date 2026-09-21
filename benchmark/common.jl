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
