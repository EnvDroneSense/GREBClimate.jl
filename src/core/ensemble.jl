"""
    run_ensemble([reduce,] run::RunSpec, configs; fields, jld2_dir="", ntasks=Threads.nthreads(),
                 logger=NullLogger(), kwargs...)

Run every [`Config`](@ref) of `configs` with [`greb_model!`](@ref), several at
a time, and return one entry per member in the order of `configs`. Each entry
is `reduce(result)`, computed in the member's task, so only what `reduce`
keeps stays in memory; without `reduce` it is the whole result.

Every member gives exactly the result of running it alone.

```julia
fields = load_climatology(dir)
configs = [preset(:co2_double; hydrology=(evaporation=e,)) for e in (:original, :skin, :skin_gust)]
warming = run_ensemble(RunSpec(ctrl=1, scnr=30), configs; fields, jld2_dir=dir) do result
    sum(global_mean(rec.Ts) for rec in result.scnr[end-11:end]) / 12
end
```

| Keyword | Meaning |
|:--------|:--------|
| `fields` | The loaded [`ClimateFields`](@ref). It is not changed: the members run on `ntasks` copies of it, each reused by one member after another |
| `jld2_dir` | The dataset directory, as for `greb_model!` |
| `ntasks` | Members running at the same time; not more than there are members. Each needs a copy of `fields`, about 370 MB |
| `logger` | Where the members' progress lines go. The default discards them; pass `Base.current_logger()` to see them |
| `kwargs` | Passed on to every `greb_model!` call, for example `allow_uninitialized`. An `observer` would be shared by all members and called from several tasks at once |

Start Julia with several threads (`-t 8,0`) for the members to run in
parallel; the gain depends on the machine. Things to keep in mind:

- A member whose physics differs from the stored flux corrections needs its
  own [`SpinUp`](@ref): the stored ones fit the configuration they were
  computed for.
- A member that diverges without going non-finite is not flagged; check the
  range in `reduce`.
- If a member throws, the other members finish and the error is rethrown.
"""
function run_ensemble(reduce, run::RunSpec, configs; fields::ClimateFields, jld2_dir::AbstractString="",
                      ntasks::Integer=Threads.nthreads(), logger=NullLogger(), kwargs...)
    ntasks >= 1 || throw(ArgumentError("ntasks must be at least 1, got $ntasks"))
    haskey(kwargs, :fields) && throw(ArgumentError("fields is given once, for all members"))
    members = collect(configs)
    isempty(members) && return []
    ntasks = min(ntasks, length(members))

    # One copy of the fields per running member. A run restores the arrays it
    # overwrites, so a copy serves one member after another.
    pool = Channel{ClimateFields}(ntasks)
    for _ in 1:ntasks
        put!(pool, deepcopy(fields))
    end

    results = Vector{Any}(undef, length(members))
    try
        @sync for (i, config) in enumerate(members)
            Threads.@spawn begin
                member_fields = take!(pool)
                try
                    results[i] = with_logger(logger) do
                        reduce(greb_model!(run, config; jld2_dir, fields=member_fields, kwargs...))
                    end
                finally
                    put!(pool, member_fields)
                end
            end
        end
    finally
        # The finished tasks stay reachable for a while and hold the pool;
        # emptying it lets the copies be collected when the call returns.
        while isready(pool)
            take!(pool)
        end
        close(pool)
    end
    return map(identity, results)   # narrows the element type
end

run_ensemble(run::RunSpec, configs; kwargs...) = run_ensemble(identity, run, configs; kwargs...)
