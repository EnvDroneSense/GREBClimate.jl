# API Reference

The configuration types ([`Config`](@ref), [`preset`](@ref), the scenario
parts, [`resolve`](@ref)) are documented on the [Configuration](@ref) page.

## Running the model

```@autodocs
Modules = [GREBClimate]
Pages = ["core/model.jl", "core/ensemble.jl", "core/output.jl", "core/postprocess.jl", "core/budgets.jl"]
```

## Data

```@autodocs
Modules = [GREBClimate]
Pages = ["src/data.jl", "src/io.jl"]
```

## State

```@autodocs
Modules = [GREBClimate]
Pages = ["core/state.jl", "core/constants.jl", "core/calendar.jl"]
```

## Physics and forcing

```@autodocs
Modules = [GREBClimate]
Pages = ["core/tendencies.jl", "forcing/forcing.jl", "physics/radiation.jl", "physics/hydrology.jl",
         "physics/ocean.jl", "physics/circulation.jl"]
```
