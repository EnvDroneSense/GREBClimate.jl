# Registry: how visualizations, variables and regions are declared.
#
# The point of this file is that adding a visualization should never mean
# editing the notebook. The notebook renders whatever is registered here, so a
# new plot is a new file in `viz/plots/` and nothing else.

"""
    VizSpec

One visualization. `compute` turns model output into numbers; `plot` turns
those numbers into a picture. They are kept separate so you can restyle
without touching the maths, or change the maths without touching the styling.

- `kind` tells the notebook which controls are relevant
  (`:series` needs a region, `:map` needs a month, and so on).
- `compute(runs, opts) -> Any` where `runs::Vector{Pair{String,Vector}}`.
- `plot(computed, opts) -> Plots.Plot`.
"""
struct VizSpec
    key::Symbol
    label::String
    kind::Symbol
    compute::Function
    plot::Function
    help::String
end

const REGISTRY = VizSpec[]

"""
    register!(; key, label, kind, compute, plot, help="")

Add a visualization. Re-registering the same `key` replaces it, so a notebook
can be re-run without accumulating duplicates.
"""
function register!(; key::Symbol, label::AbstractString, kind::Symbol,
                   compute::Function, plot::Function, help::AbstractString="")
    spec = VizSpec(key, label, kind, compute, plot, help)
    i = findfirst(s -> s.key === key, REGISTRY)
    i === nothing ? push!(REGISTRY, spec) : (REGISTRY[i] = spec)
    spec
end

"All registered visualizations, in registration order."
vizlist() = REGISTRY

"Look up one visualization by key."
function viz(key::Symbol)
    i = findfirst(s -> s.key === key, REGISTRY)
    i === nothing && error("no visualization :$key (have: $(join([s.key for s in REGISTRY], ", ")))")
    REGISTRY[i]
end

"`key => label` pairs, ready to hand to a PlutoUI `Select`."
viz_options() = [s.key => s.label for s in REGISTRY]

# ── Variables ──────────────────────────────────────────────────────────────
# Derived from the records themselves, so a new model output field becomes
# plottable everywhere without touching this file.

"Display label and unit for known fields. Unknown fields fall back gracefully."
const VAR_INFO = Dict{Symbol,NamedTuple{(:label, :unit),Tuple{String,String}}}(
    :Ts     => (label="Surface temperature",    unit="K"),
    :Ta     => (label="Air temperature",        unit="K"),
    :To     => (label="Deep-ocean temperature", unit="K"),
    :q      => (label="Specific humidity",      unit="kg/kg"),
    :albedo => (label="Albedo",                 unit="-"),
    :ice    => (label="Ice cover",              unit="fraction"),
    :precip => (label="Precipitation",          unit="mm/day"),
    :evap   => (label="Evaporation",            unit="mm/day"),
    :qcrcl  => (label="Moisture circulation",   unit="kg/kg"),
    :sw     => (label="Net shortwave",          unit="W/m²"),
    :lw     => (label="Surface longwave",       unit="W/m²"),
    :qlat   => (label="Latent heat flux",       unit="W/m²"),
    :qsens  => (label="Sensible heat flux",     unit="W/m²"),
)

"""
    variables(records)

Every field carried by the records. Derived from the data rather than a fixed
list, so a new model output field becomes plottable everywhere for free - and
so this tool does not need to import the model package at all.
"""
variables(records) = collect(fieldnames(typeof(first(records))))

"Label/unit for `v`, falling back to the bare symbol if it is not in `VAR_INFO`."
varinfo(v::Symbol) = get(VAR_INFO, v, (label=String(v), unit=""))

"`key => \"Label (unit)\"` pairs for a `Select`."
var_options(records) = [v => (varinfo(v).unit == "" ? varinfo(v).label :
                        "$(varinfo(v).label) ($(varinfo(v).unit))") for v in variables(records)]

# ── Regions ────────────────────────────────────────────────────────────────
# A region is a function (fields) -> Matrix{Bool}. `fields` may be `nothing`
# for anything that only needs latitude.

const REGIONS = Pair{Symbol,NamedTuple{(:label, :mask),Tuple{String,Function}}}[]

"""
    register_region!(key, label, mask)

`mask(nx, ny, fields) -> Matrix{Bool}`. Adding a region is one call. `fields`
is only for masks that need model data (land/ocean read `z_topo`); latitude
bands ignore it.
"""
function register_region!(key::Symbol, label::AbstractString, mask::Function)
    i = findfirst(p -> p.first === key, REGIONS)
    entry = key => (label=String(label), mask=mask)
    i === nothing ? push!(REGIONS, entry) : (REGIONS[i] = entry)
    entry
end

function region_mask(key::Symbol, nx::Int, ny::Int; fields=nothing)
    i = findfirst(p -> p.first === key, REGIONS)
    i === nothing && error("no region :$key (have: $(join([p.first for p in REGIONS], ", ")))")
    REGIONS[i].second.mask(nx, ny, fields)
end

region_options() = [p.first => p.second.label for p in REGIONS]
region_keys() = [p.first for p in REGIONS]
