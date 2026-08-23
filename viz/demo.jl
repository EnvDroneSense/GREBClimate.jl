# Runs GREB and renders every registered visualization to PNG.
#
#   julia --project=viz viz/demo.jl [output_dir]
#
# `viz/` has its own environment (Plots + a dev'd GREBClimate), so the package
# itself gains no dependency from any of this.

using Plots
gr()

using GREBClimate

include(joinpath(@__DIR__, "GREBViz.jl"))
using .GREBViz

outdir = length(ARGS) >= 1 ? ARGS[1] : pwd()
mkpath(outdir)

dir = GREBClimate.greb_data_dir(allow_download = false)
fields = load_greb_jld2!(dir; dataset = :ncep)
cfg = create_experiment_config(:co2_double)

@info "running GREB (5 yr control, 15 yr 2xCO2 scenario)"
res = redirect_stdout(devnull) do
    greb_model!(RunSpec(flux = 0, ctrl = 5, scnr = 15), cfg; jld2_dir = dir, fields = fields)
end

# `greb_runs` tags each run: control is :absolute, scenario is :anomaly.
# Getting this wrong is the difference between a readable plot and a flat line.
runs = greb_runs(res; scnr = "2xCO2")
println("runs: ", join(["$(r.label) [:$(r.kind), $(length(r.records)) months]" for r in runs], "  "))

save(p, name) = (savefig(p, joinpath(outdir, "greb_$name.png")); println("  -> greb_$name.png"))

# Mixed kinds: series-type plots split into an absolute panel and an anomaly one.
for spec in vizlist()
    spec.key === :diffmap && continue    # needs matching kinds; handled below
    save(render(spec.key, runs; var = :Ts, region = :global, month = :mean, fields = fields),
         String(spec.key))
end

# Reconstructing absolutes lets the scenario be compared with its control
# directly - in real Kelvin, and as a difference map showing WHERE it warms.
abs2x = to_absolute(runs[2], runs[1])
save(render(:timeseries, [runs[1], abs2x]; var = :Ts, annual = true), "timeseries_absolute")
save(render(:diffmap,    [abs2x, runs[1]]; var = :Ts, month = :mean, fields = fields), "diffmap_absolute")

# The reductions work without plotting anything.
println("\nfinal-year 2xCO2 Ts anomaly by region:")
for r in [:global, :tropics, :arctic, :antarctic, :land, :ocean]
    s = series(res.scnr, :Ts; region = r, fields = fields)
    println("  ", rpad(String(r), 11), round(sum(s[end-11:end]) / 12, digits = 3), " K")
end
