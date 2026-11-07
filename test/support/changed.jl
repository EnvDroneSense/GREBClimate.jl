# Run the tests a change needs, from `git diff` (staged, unstaged and untracked
# files, against HEAD). Without --run it only prints the files.
#
#   julia --project=. test/support/changed.jl [--run] [--base=REF]
#
# The map is `SOURCE_MAP` in support/selection.jl. A path it does not name runs the
# full suite. A subset proves only the files it ran: report it that way, and run
# the full suite before committing.

include(joinpath(@__DIR__, "selection.jl"))

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function changed_paths(base)
    git(args...) = readlines(Cmd(`git $args`; dir=ROOT))
    return unique(vcat(git("diff", "--name-only", base), git("ls-files", "--others", "--exclude-standard")))
end

function main(args)
    run_it = "--run" in args
    base = "HEAD"
    for a in args
        m = match(r"^--base=(.+)$", a)
        m === nothing || (base = m.captures[1])
        (a == "--run" || m !== nothing) || error("unknown argument `$a`; use --run or --base=REF")
    end
    paths = changed_paths(base)
    files, why = tests_for(paths)
    println(length(paths), " changed paths -> ", why)
    println(isempty(files) ? "no test needed" : join(files, " "))
    (run_it && !isempty(files)) || return
    code = "using Pkg; Pkg.test(test_args = $(repr(files)))"
    run(`$(Base.julia_cmd()) --project=$ROOT -e $code`)
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main(ARGS)
end
