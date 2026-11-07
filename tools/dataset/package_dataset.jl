# Build the distributable dataset packs and their SHA256. MAINTAINER TOOL.
#
# Writes one `greb_input_data-v<DATASET_VERSION>-<pack>.tar.xz` per pack (core,
# era, cmip5, enso, solar) and prints the `DATA_PACKS` block to paste into
# `src/data.jl`. Needs GNU `tar` and `xz`.
#
#   julia --project=. tools/dataset/package_dataset.jl [dataset_dir] [output_dir]
#   # defaults: greb_input_data, the current directory
#
# The archive depends only on the contents: entries are sorted, owner and group
# are zeroed, the timestamp is pinned and xz runs on one thread. Another xz
# version can still compress differently, so the checksum that counts is the one
# of the file you upload; never rebuild an archive after publishing its checksum.

using JLD2
using SHA

include(joinpath(@__DIR__, "fields.jl"))

# Stamped on every archive entry. Any constant works, but it must never change,
# or the recorded checksums stop reproducing.
const ARCHIVE_MTIME = "2020-01-01 00:00:00 UTC"

archive_name(pack::Symbol) = "greb_input_data-v$DATASET_VERSION-$pack.tar.xz"

"All `.jld2` files under `dir` as sorted paths relative to `dir`, with `/` as separator."
dataset_files(dir::AbstractString) = sort!([replace(relpath(joinpath(root, f), dir), '\\' => '/')
                                            for (root, _, names) in walkdir(dir) for f in names if endswith(f, ".jld2")])

# Keys every file carries, and the extra keys of a single-field file
const COMMON_KEYS = ("dataset_version", "description", "units", "source")
const FIELD_KEYS = ("data", "dim_names", "coords", "period")

"""
Problems with the metadata of the file at `path` (`name` is its name without ending):
missing keys, a wrong `dataset_version`, axis values that do not match the array.
"""
function metadata_problems(path::AbstractString, name::AbstractString)
    problems = String[]
    jldopen(path, "r") do file
        for key in (name in COMBINED_FILE_NAMES ? COMMON_KEYS : (COMMON_KEYS..., FIELD_KEYS...))
            haskey(file, key) || push!(problems, "lacks \"$key\"")
        end
        haskey(file, "dataset_version") && file["dataset_version"] != DATASET_VERSION &&
            push!(problems, "has dataset_version $(file["dataset_version"]), expected $DATASET_VERSION")
        if haskey(file, "data") && haskey(file, "coords") && file["data"] isa AbstractArray
            for (axis, values) in file["coords"]
                size(file["data"], axis) == length(values) ||
                    push!(problems, "axis $axis has $(length(values)) values for $(size(file["data"], axis)) points")
            end
        end
    end
    return problems
end

"""
Check `dir` against the file list: every field the model reads and every combined
file must be there, nothing else, and each file must carry its metadata.
Returns the file count and the total bytes.
"""
function validate_dataset(dir::AbstractString)
    isdir(dir) || error("dataset directory not found: $dir")
    files = dataset_files(dir)
    names = Dict(first(splitext(basename(f))) => f for f in files)
    expected = union(MODEL_FIELD_NAMES, COMBINED_FILE_NAMES)

    absent = setdiff(expected, keys(names))
    isempty(absent) || error("dataset is missing $(length(absent)) file(s):\n  ", join(sort(collect(absent)), "\n  "))
    extra = setdiff(keys(names), expected)
    isempty(extra) || error("dataset holds $(length(extra)) file(s) the model never reads:\n  ",
        join(sort(collect(extra)), "\n  "),
        "\nRegenerate without --all, or add them to MODEL_FIELD_NAMES if they are now used.")

    bad = ["$f: $p" for (name, f) in names for p in metadata_problems(joinpath(dir, f), name)]
    isempty(bad) || error("metadata problems:\n  ", join(sort(bad), "\n  "))
    return length(files), sum(f -> stat(joinpath(dir, f)).size, files)
end

"Build the archive `out` of the files `files` (relative to `dir`): a reproducible tar piped to `xz -6 -T1`."
function build_archive(dir::AbstractString, files::Vector{String}, out::AbstractString)
    list = tempname()
    try
        write(list, join(files, "\n") * "\n")
        tar = `tar --sort=name --owner=0 --group=0 --numeric-owner --mtime=$ARCHIVE_MTIME -cf - -C $dir -T $list`
        open(io -> run(pipeline(pipeline(tar, `xz -6 -T1 -c`); stdout=io)), out, "w")
    finally
        rm(list; force=true)
    end
    return out
end

function main(dir::AbstractString, outdir::AbstractString)
    nfiles, nbytes = validate_dataset(dir)
    println("$dir: $nfiles .jld2 files, $(round(nbytes / 1048576; digits=1)) MB unpacked")
    files = dataset_files(dir)
    mkpath(outdir)

    block = ["const DATA_PACKS = ("]
    sums = String[]
    for pack in PACK_NAMES
        packed = filter(f -> pack_of(f) === pack, files)
        out = build_archive(dir, packed, joinpath(outdir, archive_name(pack)))
        bytes, hash = stat(out).size, bytes2hex(open(sha256, out))
        println(rpad(basename(out), 40), lpad(bytes, 11), " bytes  ", lpad(length(packed), 2), " files  ", hash)
        push!(sums, "$hash  $(basename(out))\n")
        push!(block, "    $pack=(archive=\"$(basename(out))\", bytes=$bytes,", "        sha256=\"$hash\"),")
    end
    write(joinpath(outdir, "SHA256SUMS"), join(sums))

    println("\nPaste into src/data.jl (replacing DATA_PACKS), update the sizes in the pack messages,")
    println("and attach the archives to the GitHub release under exactly these names:\n")
    println(join(block, "\n"), "\n)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "..", "greb_input_data"),
        length(ARGS) >= 2 ? ARGS[2] : pwd())
end
