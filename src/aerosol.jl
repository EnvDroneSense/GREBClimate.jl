# Stratospheric aerosol add-on: eruptions, sustained injections and published
# records -> optical depth -> per-latitude shortwave transmission. Included
# before core/config.jl, whose PhysicsConfig holds an AerosolScenario.

"Injection classes: where the aerosol is placed latitudinally."
const AEROSOL_CLASSES = (:tropical, :tropical_nh, :tropical_sh, :nh_extratropical, :sh_extratropical)
const _TROPICAL_CLASSES = (:tropical, :tropical_nh, :tropical_sh)

# Global-mean optical depth per Tg S, back-derived from Crowley and Unterman
# (2013) and Pinatubo's 9 Tg S (Toohey and Sigl 2017); plausible 0.015-0.020.
const AOD_PER_TG_S = 0.015

# 30 / 82.1: calibrates to -30 W/m2 per unit optical depth
# (tools/diagnostics/aerosol_forcing_per_aod.jl).
const AEROSOL_TAU_SCALE = 0.365

# Latitude shapes as (latitude, value) nodes at band centres, linear between,
# constant beyond; checked by tools/diagnostics/aerosol_profile_check.jl.
# Tropical: fitted to Pinatubo's 36-month band totals in Sato-Lacis and GloSSAC.
const _TROPICAL_NODES = ((7.5, 1.0), (22.5, 0.68), (37.5, 0.65), (52.5, 0.65), (75.0, 0.55))
# One-hemisphere tropical, latitude positive in the eruption's hemisphere:
# fitted to El Chichon and Agung the same way.
const _ONE_HEMISPHERE_NODES = ((-22.5, 0.28), (-7.5, 0.62), (7.5, 1.0), (22.5, 1.0), (37.5, 0.75))
# Extratropical, eruption's hemisphere only: Katmai (Stothers 1996), nothing
# south of 30 N, 0.13 / 0.23 at 30-45 N; polar value assumed equal to 45-60 N.
const _EXTRATROPICAL_NODES = ((30.0, 0.0), (37.5, 0.57), (52.5, 1.0))

# Extratropical (rise, decay) in years: Katmai's 0.8-year decay and a peak at
# 2.5 months (Stothers 1996). The tropical classes use `tau_rise`/`tau_decay`.
const _EXTRATROPICAL_TAU = (rise = 0.082, decay = 0.8)

function _interp(nodes, x::Real)
    x <= first(nodes)[1] && return first(nodes)[2]
    for k in 2:length(nodes)
        (x1, v1), (x0, v0) = nodes[k], nodes[k - 1]
        x < x1 && return v0 + (v1 - v0) * (x - x0) / (x1 - x0)
    end
    return last(nodes)[2]
end

_tropical_shape(phi::Real) = _interp(_TROPICAL_NODES, abs(phi))
_extratropical_shape(a::Real) = _interp(_EXTRATROPICAL_NODES, a)

_check_class(cls::Symbol) = cls in AEROSOL_CLASSES ||
    throw(ArgumentError("unknown injection class :$cls; valid: $AEROSOL_CLASSES"))

function _class_profile(cls::Symbol)
    h = cls in (:tropical_nh, :nh_extratropical) ? 1.0 : -1.0
    w = if cls === :tropical
        [_tropical_shape(Float64(phi)) for phi in lat_grid]
    elseif cls in (:tropical_nh, :tropical_sh)
        [_interp(_ONE_HEMISPHERE_NODES, h * Float64(phi)) for phi in lat_grid]
    else
        [_extratropical_shape(h * Float64(phi)) for phi in lat_grid]
    end
    c = [cosd(Float64(phi)) for phi in lat_grid]
    return w ./ (sum(w .* c) / sum(c))
end

# Latitude profile per class, scaled to an area-weighted mean of 1. Read-only.
const _CLASS_PROFILES = Dict{Symbol,Vector{Float64}}(cls => _class_profile(cls) for cls in AEROSOL_CLASSES)

"""
    Eruption(year, cls; tg_s=nothing, peak_aod=nothing, aod_per_tg_s=0.015)

One volcanic injection at decimal `year` on the model clock (the scenario run
starts at 1950 for most experiments). `cls` sets where the aerosol goes:

- `:tropical`: both hemispheres, like Pinatubo;
- `:tropical_nh`, `:tropical_sh`: mostly one hemisphere, like El Chichon, Agung;
- `:nh_extratropical`, `:sh_extratropical`: poleward of 30 degrees in one
  hemisphere, like Katmai (eruptions beyond 25 degrees latitude).

Give exactly one of `tg_s` (injected sulfur in Tg S, times `aod_per_tg_s`) or
`peak_aod` (global-mean peak optical depth at 550 nm). The linear mass scaling
overestimates eruptions much larger than Pinatubo.
"""
struct Eruption
    year::Float64
    cls::Symbol
    peak_aod::Float64
    function Eruption(year::Real, cls::Symbol, peak_aod::Real)
        _check_class(cls)
        peak_aod >= 0 || throw(ArgumentError("peak_aod must be non-negative, got $peak_aod"))
        return new(Float64(year), cls, Float64(peak_aod))
    end
end

function Eruption(year::Real, cls::Symbol; tg_s::Union{Real,Nothing}=nothing,
                  peak_aod::Union{Real,Nothing}=nothing, aod_per_tg_s::Real=AOD_PER_TG_S)
    (tg_s === nothing) != (peak_aod === nothing) ||
        throw(ArgumentError("give exactly one of tg_s or peak_aod"))
    return Eruption(year, cls, peak_aod === nothing ? tg_s * aod_per_tg_s : peak_aod)
end

"""
    SustainedInjection(start_year, stop_year, cls, target_aod)

Holds the global-mean optical depth at `target_aod` between `start_year` and
`stop_year` (use `Inf` for no stop). The optical depth rises and, after the
stop, decays with the timescales of its class (see [`AerosolScenario`](@ref));
choosing a stop year inside the run gives a termination shock.
"""
struct SustainedInjection
    start_year::Float64
    stop_year::Float64
    cls::Symbol
    target_aod::Float64
    function SustainedInjection(start_year::Real, stop_year::Real, cls::Symbol, target_aod::Real)
        _check_class(cls)
        stop_year > start_year ||
            throw(ArgumentError("stop_year must exceed start_year, got $start_year and $stop_year"))
        target_aod >= 0 || throw(ArgumentError("target_aod must be non-negative, got $target_aod"))
        return new(Float64(start_year), Float64(stop_year), cls, Float64(target_aod))
    end
end

"""
    AerosolSeries

A published optical-depth record already interpolated onto the model latitudes:
`years` (ascending decimal years on the model calendar), `aod` of size
`(ydim, length(years))` so one time's column is contiguous, and `source`, the
provenance string stored in the file. Build one with [`load_aerosol_series`](@ref).
"""
struct AerosolSeries
    years::Vector{Float64}
    aod::Matrix{Float32}
    source::String
end

# `load_aerosol_series` warns when a record ends above this optical depth.
const SERIES_END_WARN_AOD = 0.005

"""
    AerosolScenario(; eruptions=[], injections=[], series=nothing, tau_rise=0.21,
                    tau_decay=1.0, ssa=1.0, asymmetry=0.7, mu0=0.5,
                    tau_scale=0.365)

Stratospheric aerosol add-on for `cfg.aerosol`: dims the shortwave per
latitude in the scenario run; the control run stays aerosol-free.

Optical depth is the sum of `eruptions`, `injections` and `series` (a published
record, see [`load_aerosol_series`](@ref); an eruption also in the record counts
twice). Eruptions and injections rise and decay through two linear reservoirs:
for the tropical classes with `tau_rise` and `tau_decay` (years; the defaults
peak at 5 months and decay in a year, like Pinatubo), for the extratropical
classes with a 2.5-month peak and a 0.8-year decay (Katmai).

Transmission is delta-Eddington with `ssa`, `asymmetry` and `mu0` on `tau_scale`
times the optical depth. The default `tau_scale` gives about -30 W/m2 per unit
optical depth (Sato et al. 1993) for the default `ssa`, `asymmetry` and `mu0`;
if you change those, recalibrate it with
`tools/diagnostics/aerosol_forcing_per_aod.jl`. Longwave effects, stratospheric
heating, ozone chemistry and particle growth are not represented.
"""
struct AerosolScenario
    eruptions::Vector{Eruption}
    injections::Vector{SustainedInjection}
    series::Union{Nothing,AerosolSeries}
    tau_rise::Float64
    tau_decay::Float64
    ssa::Float64
    asymmetry::Float64
    mu0::Float64
    tau_scale::Float64
    # Per class: (rise, decay, 1 / pulse peak), derived from the fields above.
    chains::NamedTuple{AEROSOL_CLASSES,NTuple{length(AEROSOL_CLASSES),NTuple{3,Float64}}}
end

function AerosolScenario(eruptions, injections, series, tau_rise, tau_decay, ssa, asymmetry, mu0, tau_scale)
    chain(tr, td) = (tr, td, 1 / chain_impulse(chain_peak_time(tr, td), tr, td))
    chains = map(AEROSOL_CLASSES) do cls
        cls in _TROPICAL_CLASSES ? chain(tau_rise, tau_decay) : chain(_EXTRATROPICAL_TAU.rise, _EXTRATROPICAL_TAU.decay)
    end
    return AerosolScenario(eruptions, injections, series, tau_rise, tau_decay, ssa, asymmetry, mu0,
                           tau_scale, NamedTuple{AEROSOL_CLASSES}(chains))
end

function AerosolScenario(; eruptions=Eruption[], injections=SustainedInjection[],
        series::Union{Nothing,AerosolSeries}=nothing, tau_rise::Real=0.21,
        tau_decay::Real=1.0, ssa::Real=1.0, asymmetry::Real=0.7, mu0::Real=0.5,
        tau_scale::Real=AEROSOL_TAU_SCALE)
    0 < tau_rise < tau_decay ||
        throw(ArgumentError("need 0 < tau_rise < tau_decay (years), got $tau_rise and $tau_decay"))
    0 <= ssa <= 1 || throw(ArgumentError("ssa must be in [0, 1], got $ssa"))
    0 <= asymmetry < 1 || throw(ArgumentError("asymmetry must be in [0, 1), got $asymmetry"))
    0 < mu0 <= 1 || throw(ArgumentError("mu0 must be in (0, 1], got $mu0"))
    0 <= tau_scale < Inf || throw(ArgumentError("tau_scale must be finite and >= 0, got $tau_scale"))
    return AerosolScenario(collect(Eruption, eruptions), collect(SustainedInjection, injections),
        series, Float64(tau_rise), Float64(tau_decay), Float64(ssa), Float64(asymmetry), Float64(mu0),
        Float64(tau_scale))
end

# ── Two-reservoir chain (times in years) ─────────────────────────────────

"Response of the aerosol reservoir to a unit mass injected at t = 0."
function chain_impulse(t, tr, td)
    t <= 0 && return 0.0
    return td / (td - tr) * (exp(-t / td) - exp(-t / tr))
end

"Time of the maximum of [`chain_impulse`](@ref)."
chain_peak_time(tr, td) = log(td / tr) * tr * td / (td - tr)

"Response to a constant source switched on at t = 0, scaled to reach 1."
function chain_step(t, tr, td)
    t <= 0 && return 0.0
    return 1.0 - (td * exp(-t / td) - tr * exp(-t / tr)) / (td - tr)
end

"""
    aerosol_optical_depth!(aod::Vector{Float32}, sc::AerosolScenario, t) -> aod

Fills `aod` (length `ydim`) with the optical depth at decimal year `t`.
Depends only on `t`; allocates nothing.
"""
function aerosol_optical_depth!(aod::Vector{Float32}, sc::AerosolScenario, t::Real)
    fill!(aod, 0.0f0)
    for e in sc.eruptions
        tr, td, inv_peak = sc.chains[e.cls]    # inv_peak scales the pulse to a peak of 1
        a = e.peak_aod * inv_peak * chain_impulse(t - e.year, tr, td)
        a == 0.0 && continue
        aod .+= a .* _CLASS_PROFILES[e.cls]
    end
    for s in sc.injections
        tr, td, _ = sc.chains[s.cls]
        u = chain_step(t - s.start_year, tr, td)
        isfinite(s.stop_year) && (u -= chain_step(t - s.stop_year, tr, td))
        u == 0.0 && continue
        aod .+= (s.target_aod * u) .* _CLASS_PROFILES[s.cls]
    end
    sc.series === nothing || _add_series!(aod, sc.series, t)
    return aod
end

# Adds the record at decimal year `t`, linear in time, zero outside its range.
function _add_series!(aod::Vector{Float32}, s::AerosolSeries, t::Real)
    yrs = s.years
    (t < first(yrs) || t > last(yrs)) && return aod
    k = searchsortedlast(yrs, t)
    if k == length(yrs)
        aod .+= @view s.aod[:, k]
    else
        frac = (t - yrs[k]) / (yrs[k + 1] - yrs[k])
        aod .+= (1 - frac) .* @view(s.aod[:, k]) .+ frac .* @view(s.aod[:, k + 1])
    end
    return aod
end

# ── Delta-Eddington transmission ─────────────────────────────────────────
# Joseph, Wiscombe and Weinman (1976): the Eddington two-stream solution with
# delta-scaled parameters. Closed forms: Meador and Weaver (1980), eqs. 14-18
# and Table 1.

"""
    delta_eddington(tau, ssa, g, mu0) -> (R, T)

Reflectance `R` and total (direct plus diffuse) transmission `T` of a
scattering layer over a black surface, by the delta-Eddington two-stream
method. `tau` is the optical depth, `ssa` the single-scattering albedo, `g`
the asymmetry parameter and `mu0` the cosine of the solar zenith angle. Computed
in `Float64`: conservative scattering (`ssa = 1`) makes the closed form
singular, so the scaled albedo is clamped just below 1.
"""
function delta_eddington(tau::Real, ssa::Real, g::Real, mu0::Real)
    tau <= 0 && return (0.0, 1.0)
    tau / mu0 > 300 &&
        throw(DomainError(tau, "tau / mu0 above 300 overflows the two-stream exponentials"))
    f = g^2
    τ = (1 - ssa * f) * tau
    ω = min((1 - f) * ssa / (1 - ssa * f), 1 - 1e-9)
    gs = g / (1 + g)
    γ1 = (7 - ω * (4 + 3gs)) / 4
    γ2 = -(1 - ω * (4 - 3gs)) / 4
    k = sqrt(γ1^2 - γ2^2)
    # k * mu0 = 1 zeroes D below; the singularity is removable, so step off it.
    mu = abs(1 - (k * mu0)^2) < 1e-8 ? mu0 * (1 + 1e-6) : mu0
    γ3 = (2 - 3gs * mu) / 4
    γ4 = 1 - γ3
    α1 = γ1 * γ4 + γ2 * γ3
    α2 = γ1 * γ3 + γ2 * γ4
    ekt, emkt, emt = exp(k * τ), exp(-k * τ), exp(-τ / mu)
    D = (1 - (k * mu)^2) * ((k + γ1) * ekt + (k - γ1) * emkt)
    R = ω / D * ((1 - k * mu) * (α2 + k * γ3) * ekt -
                 (1 + k * mu) * (α2 - k * γ3) * emkt -
                 2k * (γ3 - α2 * mu) * emt)
    T = emt * (1 - ω / D * ((1 + k * mu) * (α1 + k * γ4) * ekt -
                            (1 - k * mu) * (α1 - k * γ4) * emkt -
                            2k * (γ4 + α1 * mu) / emt))  # / emt is exp(+τ/μ0) in eq. 15
    return (R, T)
end

"""
    aerosol_transmission!(mult::Vector{Float32}, aod::Vector{Float32}, sc::AerosolScenario, t) -> mult

Multiplies `mult` (length `ydim`) by the aerosol shortwave transmission at
decimal year `t`, using `aod` as scratch. Rows without aerosol are untouched, so
the multiplier is bit-identical when no aerosol is present.
"""
function aerosol_transmission!(mult::Vector{Float32}, aod::Vector{Float32}, sc::AerosolScenario, t::Real)
    aerosol_optical_depth!(aod, sc, t)
    for j in 1:ydim
        tau = sc.tau_scale * Float64(aod[j])
        tau == 0.0 && continue
        mult[j] *= Float32(delta_eddington(tau, sc.ssa, sc.asymmetry, sc.mu0)[2])
    end
    return mult
end
