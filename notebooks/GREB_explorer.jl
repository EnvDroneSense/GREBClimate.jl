### A Pluto.jl notebook ###
# v1.0.3

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

# ╔═╡ b37f2536-200b-49b1-a823-d5625431c1c7
md"""
# GREB explorer

Run an experiment, then look at the output with the `viz/` plots.
"""

# ╔═╡ 2759e19a-df8b-45d3-9d66-f45b9498f972
md"""## 1. Run"""

# ╔═╡ 34ed8352-d79c-418d-ae5b-014b09a8ff30
begin
    data_dir = greb_data_dir(allow_download=false)
    fields = load_greb_jld2!(data_dir; dataset=:ncep)
end;

# ╔═╡ 81e2898f-92b5-4bfa-837d-82319d446a7b
md"""Experiment $(@bind experiment Select(sort(collect(keys(GREBClimate._EXPERIMENT_OVERRIDES))), default=:co2_double))"""

# ╔═╡ 5e5b94a6-59d6-4f92-9978-c0ca505a5f2e
md"""Control years $(@bind n_ctrl Slider(1:30, default=5, show_value=true))

Scenario years $(@bind n_scnr Slider(0:50, default=15, show_value=true))

Run the model $(@bind go CheckBox(default=false))"""

# ╔═╡ 79fdb975-0f02-43f7-9a05-39e8f5408c5f
result = go ? redirect_stdout(devnull) do
    greb_model!(RunSpec(flux=3, ctrl=n_ctrl, scnr=n_scnr), create_experiment_config(experiment);
                jld2_dir=data_dir, fields=deepcopy(fields))   # the model mutates fields
end : nothing;

# ╔═╡ 372f4e9a-f1be-4e57-b293-b2592e29c60f
md"""## 2. Plots"""

# ╔═╡ df41b455-c829-4ae6-bd04-263f0b28a9e4
md"""
Variable $(@bind v Select([k => fieldinfo(k).label for k in (result === nothing ? [:Ts] : keys(first(result.ctrl)))]))

Map month $(@bind mon Select([:mean => "run mean", :last => "final month"]))

Annual means $(@bind yearly CheckBox(default=false))
"""

# ╔═╡ fba7c888-e783-442e-8d1a-ab66a7629d8c
result === nothing ? md"*Tick **Run the model** above.*" : plot_map(result; var=v, month=mon, fields=fields)

# ╔═╡ 8b90b23a-8014-40fc-9fb0-1b3c04a30f53
result === nothing ? md"" : plot_timeseries(result; var=v, annual=yearly)

# ╔═╡ a0e60d10-c2f8-42ae-9cb5-e7b8b4787b7e
result === nothing ? md"" : plot_seasonal(result; var=v)

# ╔═╡ ca0cefe0-ac19-4e8c-ac4c-0ff3075ccbfd
result === nothing ? md"" : plot_hovmoller(result; var=v)

# ╔═╡ 604877f9-26bc-4dde-abc4-f7969831d3f7
md"""## 3. Through the run

Step $(@bind evo_step Select([:year => "yearly means", :month => "monthly"]))
"""

# ╔═╡ 12f808c0-c3c3-4154-ad1b-3abcee55658f
evo = result === nothing ? nothing : evolution(result; var=v, step=evo_step);

# ╔═╡ aff47e88-daa7-40f7-9167-5ee510419968
md"""
Frame $(@bind evo_frame Slider(1:(evo === nothing ? 1 : frame_count(evo)), show_value=true))

Play $(@bind evo_play CheckBox(default=false)) $(@bind evo_tick Clock(0.4))

*To animate, tick **Play** and press **Start**.*
"""

# ╔═╡ 65a00347-bfe3-4445-848d-9d9a80d1b4c5
if evo === nothing
    md""
else
    # The clock never resets, so it picks the frame only while Play is ticked.
    evolution_frame(evo, evo_play ? mod1(evo_tick, frame_count(evo)) : evo_frame; fields=fields)
end

# ╔═╡ 06071081-4a67-4eb6-a30d-dc10460457b3
md"""
Write GIF $(@bind evo_write CheckBox(default=false))

frames/s $(@bind evo_fps Slider(2:12, default=6, show_value=true))
"""

# ╔═╡ 90176953-2498-4b7e-94d9-89723cf3322e
if evo === nothing || !evo_write
    md"*Tick **Write GIF** to render every frame (a 15-year monthly run takes ~15 s).*"
else
    let path = joinpath(tempdir(), "greb_evolution_$(evo.var)_$(evo.step).gif")
        evolution_gif(path, evo; fps=evo_fps, fields=fields)
        md"$(DownloadButton(read(path), basename(path))) $(LocalResource(path))"
    end
end

# ╔═╡ Cell order:
# ╟─b37f2536-200b-49b1-a823-d5625431c1c7
# ╟─5f08d8d1-c39d-4bb4-bead-7133854f0715
# ╟─2759e19a-df8b-45d3-9d66-f45b9498f972
# ╟─34ed8352-d79c-418d-ae5b-014b09a8ff30
# ╟─81e2898f-92b5-4bfa-837d-82319d446a7b
# ╟─5e5b94a6-59d6-4f92-9978-c0ca505a5f2e
# ╟─79fdb975-0f02-43f7-9a05-39e8f5408c5f
# ╟─372f4e9a-f1be-4e57-b293-b2592e29c60f
# ╟─df41b455-c829-4ae6-bd04-263f0b28a9e4
# ╟─fba7c888-e783-442e-8d1a-ab66a7629d8c
# ╟─8b90b23a-8014-40fc-9fb0-1b3c04a30f53
# ╟─a0e60d10-c2f8-42ae-9cb5-e7b8b4787b7e
# ╟─ca0cefe0-ac19-4e8c-ac4c-0ff3075ccbfd
# ╟─604877f9-26bc-4dde-abc4-f7969831d3f7
# ╟─12f808c0-c3c3-4154-ad1b-3abcee55658f
# ╟─aff47e88-daa7-40f7-9167-5ee510419968
# ╟─65a00347-bfe3-4445-848d-9d9a80d1b4c5
# ╟─06071081-4a67-4eb6-a30d-dc10460457b3
# ╟─90176953-2498-4b7e-94d9-89723cf3322e
