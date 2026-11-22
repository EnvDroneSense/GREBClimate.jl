# =============================================================================
# check_docs.jl - static documentation checks for GREBClimate.jl.
#
# Dependency-free (stdlib only). Does NOT load the package, so it is fast and
# cannot be broken by a compile error:
#
#   julia .claude/skills/docs-check/check_docs.jl
#
# Checks:
#   1. Markdown links to local files resolve, relative to the linking file.
#   2. Julia version claims agree across README, docs and Project.toml.
#
# Exits nonzero if any check fails, so it gates CI. Exported names and
# docstrings are checked by the docs build (`checkdocs = :exports`).
#
# It deliberately does NOT run the documented snippets - that needs the real
# package and dataset. See SKILL.md; running them is the check that matters
# most and has no substitute.
# =============================================================================

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
# Folders that hold no pages of this repository: other environments, build output, worktrees.
const SKIP_DIRS = Set([".git", ".julia", ".claude/worktrees", "viz", "docs/build", "dataset_packs",
    "greb_input_data", "Data", "benchmark/profiles", "node_modules"])

function markdown_files()
    files = String[]
    for (dir, dirs, names) in walkdir(ROOT)
        rel = replace(relpath(dir, ROOT), '\\' => '/')
        filter!(d -> !((rel == "." ? d : "$rel/$d") in SKIP_DIRS), dirs)
        for n in names
            endswith(n, ".md") && push!(files, replace(relpath(joinpath(dir, n), ROOT), '\\' => '/'))
        end
    end
    return sort(files)
end

failures = String[]
fail(msg) = (push!(failures, msg); println("  FAIL  ", msg))
ok(msg) = println("  ok    ", msg)

# -- 1. Markdown links to local files -----------------------------------------
# Only inline links `[text](target)`; bare filenames in prose and directory
# trees are not links and are not checked.
println("\n[1] local markdown links resolve")
let checked = 0, files = markdown_files()
    for rel in files
        path = joinpath(ROOT, rel)
        base = dirname(path)
        for m in eachmatch(r"\[[^\]]*\]\(([^)\s]+)\)", read(path, String))
            target = m.captures[1]
            # skip URLs, anchors, and Documenter @ref links
            (startswith(target, "http") || startswith(target, "#") ||
             startswith(target, "mailto:") || occursin("@ref", target)) && continue
            file = first(split(target, '#'))           # drop any #anchor
            isempty(file) && continue
            checked += 1
            ispath(normpath(joinpath(base, file))) ||
                fail("$rel -> $target (does not exist)")
        end
    end
    ok("$checked local links in $(length(files)) files checked")
end

# -- 2. Julia version claims agree ---------------------------------------------
println("\n[2] Julia version claims agree")
let claims = Dict{String,Set{String}}()
    proj = read(joinpath(ROOT, "Project.toml"), String)
    m = match(r"^julia\s*=\s*\"([^\"]+)\""m, proj)
    m === nothing ? fail("Project.toml has no julia compat bound") :
    (claims["Project.toml"] = Set([m.captures[1]]))
    for rel in ["README.md", "docs/src/index.md"]
        s = read(joinpath(ROOT, rel), String)
        vs = Set(m.captures[1] for m in eachmatch(r"Julia\s+\*{0,2}(1\.\d+)"i, s))
        isempty(vs) || (claims[rel] = vs)
    end
    allv = union(values(claims)...)
    if length(allv) == 1
        ok("all sources say Julia $(first(allv))")
    else
        parts = [string(k, " => ", join(sort(collect(v)), ",")) for (k, v) in claims]
        fail("disagreement: " * join(parts, "; "))
    end
end

println()
if isempty(failures)
    println("docs-check: all static checks passed.")
    println("STILL REQUIRED: run the documented snippets. See SKILL.md section 1.")
else
    println("docs-check: $(length(failures)) failure(s).")
    exit(1)
end
