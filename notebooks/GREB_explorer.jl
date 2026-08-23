### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ b37f2536-200b-49b1-a823-d5625431c1c7
md"""
# GREB explorer

Run an experiment, then explore the output. Every plot below comes from the
`viz/` registry - drop a file in `viz/plots/` and it appears in the menu with
no change to this notebook.
"""

# ╔═╡ 5f08d8d1-c39d-4bb4-bead-7133854f0715
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "viz"))
    using Plots, PlutoUI, GREBClimate
    include(joinpath(@__DIR__, "..", "viz", "GREBViz.jl"))
    using .GREBViz
    gr()
    md"*environment ready*"
end

# ╔═╡ 2759e19a-df8b-45d3-9d66-f45b9498f972
md"""## 1. Data"""

# ╔═╡ 34ed8352-d79c-418d-ae5b-014b09a8ff30
fields_master = load_greb_jld2!(GREBClimate.greb_data_dir(allow_download=false); dataset=:ncep);

# ╔═╡ 504a2ed8-f4c1-4607-85cc-8b0c68117f70
md"""## 2. Experiment"""

# ╔═╡ 81e2898f-92b5-4bfa-837d-82319d446a7b
begin
    # Prefer the package's own experiment table so this menu never drifts from
    # what the model actually supports. It is a private const, so fall back to a
    # curated list rather than breaking if it is ever renamed.
    experiment_choices = isdefined(GREBClimate, :_EXPERIMENT_OVERRIDES) ?
        sort(collect(keys(GREBClimate._EXPERIMENT_OVERRIDES))) :
        [:full_model, :co2_double, :co2_quadruple, :co2_half, :rcp85,
         :elnino, :lanina, :sst_plus1, :solar_plus27, :paleo_231kyr]
    md"""Experiment $(@bind experiment Select(experiment_choices, default=:co2_double))"""
end

# ╔═╡ 5e5b94a6-59d6-4f92-9978-c0ca505a5f2e
md"""Control years $(@bind n_ctrl Slider(1:30, default=5, show_value=true)) &nbsp;
Scenario years $(@bind n_scnr Slider(0:50, default=15, show_value=true))"""

# ╔═╡ 2589f0ec-b842-4728-82df-86273b8d0803
md"""Tick to run the model (it re-runs whenever anything above changes):
$(@bind go CheckBox(default=false))"""

# ╔═╡ 79fdb975-0f02-43f7-9a05-39e8f5408c5f
result = if go
    cfg = create_experiment_config(experiment)
    f = deepcopy(fields_master)
    redirect_stdout(devnull) do
        greb_model!(RunSpec(flux=0, ctrl=n_ctrl, scnr=n_scnr), cfg;
                    jld2_dir=GREBClimate.greb_data_dir(allow_download=false), fields=f)
    end
else
    nothing
end;

# ╔═╡ d5a2c14d-2b20-47e4-8d36-9a5188fda10a
runs = result === nothing ? nothing : greb_runs(result; scnr=String(experiment));

# ╔═╡ ca0cefe0-ac19-4e8c-ac4c-0ff3075ccbfd
if runs === nothing
    md"""!!! warning "Not run yet"
        Tick the box above to run the model."""
else
    md"""Ready: $(join(["**$(r.label)** (`:$(r.kind)`, $(length(r.records)) months)" for r in runs], " · "))"""
end

# ╔═╡ 372f4e9a-f1be-4e57-b293-b2592e29c60f
md"""## 3. Explore"""

# ╔═╡ df41b455-c829-4ae6-bd04-263f0b28a9e4
md"""
Plot $(@bind which Select(viz_options())) &nbsp;
Variable $(@bind v Select(var_options(runs === nothing ? [(; Ts=zeros(Float32,1,1))] : runs[1].records))) &nbsp;
Region $(@bind reg Select(region_options()))
"""

# ╔═╡ b120788d-c724-4d8c-b2dc-03a1a616a79b
md"""
Annual means $(@bind yearly CheckBox(default=false)) &nbsp;
Coastlines $(@bind coast CheckBox(default=true)) &nbsp;
Map month $(@bind mon Select([:mean => "run mean", :last => "final month"]))
"""

# ╔═╡ fba7c888-e783-442e-8d1a-ab66a7629d8c
if runs === nothing
    md"*run the model first*"
else
    try
        render(which, runs; var=v, region=reg, annual=yearly,
               month=mon, coastlines=coast, fields=fields_master)
    catch e
        md"""!!! danger "Cannot draw this"
            $(sprint(showerror, e))"""
    end
end

# ╔═╡ a0e60d10-c2f8-42ae-9cb5-e7b8b4787b7e
md"""## 4. Scenario in real units

A scenario from `greb_model!` is an *anomaly* against its control, which is why
it gets its own panel above. `to_absolute` adds the control climatology back so
the two can share an axis - and so a difference map means something.
"""

# ╔═╡ 8b90b23a-8014-40fc-9fb0-1b3c04a30f53
if runs === nothing || length(runs) < 2 || length(runs[1].records) < 12
    md"*needs a control of at least 12 months and a scenario*"
else
    abs_scnr = to_absolute(runs[2], runs[1])
    render(:diffmap, [abs_scnr, runs[1]]; var=v, month=mon,
           coastlines=coast, fields=fields_master)
end

# ╔═╡ df59a1d6-f019-442d-a9c4-a7b80f2a3c28
md"""### What is registered

Adding a file to `viz/plots/` extends this list, and the menu above, automatically.
"""

# ╔═╡ 8f4020b1-52f7-4132-9a9b-d5f3c7071347
begin
    _lines = String[]
    for spec in vizlist()
        push!(_lines, "- `:$(spec.key)` \u2014 **$(spec.label)**")
        isempty(spec.help) || push!(_lines, "   $(spec.help)")
    end
    Markdown.parse(join(_lines, "\n"))
end

# ╔═╡ Cell order:
# ╟─b37f2536-200b-49b1-a823-d5625431c1c7
# ╠═5f08d8d1-c39d-4bb4-bead-7133854f0715
# ╟─2759e19a-df8b-45d3-9d66-f45b9498f972
# ╠═34ed8352-d79c-418d-ae5b-014b09a8ff30
# ╟─504a2ed8-f4c1-4607-85cc-8b0c68117f70
# ╠═81e2898f-92b5-4bfa-837d-82319d446a7b
# ╟─5e5b94a6-59d6-4f92-9978-c0ca505a5f2e
# ╟─2589f0ec-b842-4728-82df-86273b8d0803
# ╠═79fdb975-0f02-43f7-9a05-39e8f5408c5f
# ╠═d5a2c14d-2b20-47e4-8d36-9a5188fda10a
# ╠═ca0cefe0-ac19-4e8c-ac4c-0ff3075ccbfd
# ╟─372f4e9a-f1be-4e57-b293-b2592e29c60f
# ╟─df41b455-c829-4ae6-bd04-263f0b28a9e4
# ╟─b120788d-c724-4d8c-b2dc-03a1a616a79b
# ╠═fba7c888-e783-442e-8d1a-ab66a7629d8c
# ╟─a0e60d10-c2f8-42ae-9cb5-e7b8b4787b7e
# ╠═8b90b23a-8014-40fc-9fb0-1b3c04a30f53
# ╟─df59a1d6-f019-442d-a9c4-a7b80f2a3c28
# ╠═8f4020b1-52f7-4132-9a9b-d5f3c7071347
