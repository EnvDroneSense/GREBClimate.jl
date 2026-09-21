# Run GREB and save every plot to PNG.
#
#   julia --project=viz viz/demo.jl [output_dir]    # default: <tempdir>/greb_viz_demo

using Plots, GREBClimate
gr()
include(joinpath(@__DIR__, "GREBViz.jl"))
using .GREBViz

outdir = mkpath(get(ARGS, 1, joinpath(tempdir(), "greb_viz_demo")))
dir = greb_data_dir(allow_download=false)
fields = load_greb_jld2!(dir; dataset=:ncep)

@info "running GREB (3 yr flux spin-up, 5 yr control, 15 yr 2xCO2 scenario)"
res = redirect_stdout(devnull) do
    greb_model!(RunSpec(flux=3, ctrl=5, scnr=15), create_experiment_config(:co2_double);
                jld2_dir=dir, fields=deepcopy(fields))
end

save(p, name) = (savefig(p, joinpath(outdir, "greb_$name.png")); println("  -> greb_$name.png"))
save(plot_map(res; fields=fields), "map")
save(plot_timeseries(res), "timeseries")
save(plot_timeseries(res; annual=true), "timeseries_annual")
save(plot_seasonal(res), "seasonal")
save(plot_hovmoller(res), "hovmoller")
evolution_gif(joinpath(outdir, "greb_evolution.gif"), evolution(res; step=:year); fields=fields)
println("  -> greb_evolution.gif")

s = series(res.scnr, :Ts)
println("\nfinal-year 2xCO2 global-mean Ts anomaly: ", round(sum(s[end-11:end]) / 12, digits=3), " K")
println("output in ", outdir)
