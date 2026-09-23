# Stratospheric aerosol add-on: eruptions and sustained injection -> optical
# depth -> per-latitude shortwave transmission. Included after constants.jl
# (needs `lat_grid`, `ydim`) and before config.jl (PhysicsConfig will hold an
# AerosolScenario).

"Injection classes: where the aerosol is placed latitudinally."
const AEROSOL_CLASSES = (:tropical, :nh_extratropical, :sh_extratropical)

# Global-mean visible optical depth per Tg S injected. Back-derived from
# Crowley and Unterman (2013) and a Pinatubo injection of 9 Tg S (Toohey and
# Sigl 2017); the plausible range is 0.015 to 0.020.
const AOD_PER_TG_S = 0.015

# Extratropical shapes: (centre latitude, gaussian half-width) in degrees.
# Provisional; no source values were found for them.
const _AEROSOL_EXTRATROPICAL_SHAPE = (
    nh_extratropical = (50.0, 30.0),
    sh_extratropical = (-50.0, 30.0),
)

# Tropical shape, from the four equal-area bands of Crowley and Unterman
# (2013): 1 in the low-latitude bands (centres 15 degrees), 0.8 in the
# high-latitude bands (centres 60 degrees), linear between, constant beyond.
function _tropical_shape(phi::Real)
    a = abs(phi)
    a <= 15 && return 1.0
    a >= 60 && return 0.8
    return 1.0 - 0.2 * (a - 15) / 45
end

_check_class(cls::Symbol) = cls in AEROSOL_CLASSES ||
    throw(ArgumentError("unknown injection class :$cls; valid: $AEROSOL_CLASSES"))

function _class_profile(cls::Symbol)
    w = if cls === :tropical
        [_tropical_shape(Float64(phi)) for phi in lat_grid]
    else
        (phi0, sigma) = _AEROSOL_EXTRATROPICAL_SHAPE[cls]
        [exp(-((Float64(phi) - phi0) / sigma)^2) for phi in lat_grid]
    end
    c = [cosd(Float64(phi)) for phi in lat_grid]
    return w ./ (sum(w .* c) / sum(c))
end

# Built once at load and never written afterwards, so concurrent reads are safe.
# The add-ons run on the main thread in any case; only `circulation!` spawns.
const _CLASS_PROFILES = Dict{Symbol,Vector{Float64}}(cls => _class_profile(cls) for cls in AEROSOL_CLASSES)

"""
    Eruption(year, cls; tg_s=nothing, peak_aod=nothing)

One volcanic injection at decimal `year` on the model clock (the scenario run
starts at 1950 for most experiments). `cls` is `:tropical`,
`:nh_extratropical` or `:sh_extratropical`. Give exactly one of `tg_s`
(injected sulfur, Tg S, converted with a linear scaling) or `peak_aod`
(global-mean peak visible optical depth, used as is). The linear mass scaling
holds for most eruptions and ignores the loss of efficiency above roughly
Pinatubo-sized peaks. The `:tropical` class spreads the aerosol into both
hemispheres, which suits Pinatubo but not hemispherically asymmetric tropical
eruptions such as Agung or El Chichon.
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
stop, decays with the scenario's `tau_rise` and `tau_decay`; choosing a stop
year inside the run gives a termination shock.
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
`years` (ascending decimal years), `aod` of size `(ydim, length(years))` so one
year's column is contiguous, and `source`, the provenance string stored in the
file. Build one with [`load_aerosol_series`](@ref).
"""
struct AerosolSeries
    years::Vector{Float64}
    aod::Matrix{Float32}
    source::String
end

# Loading warns when the last year's largest optical depth exceeds this. Chosen
# as a small fraction of a Pinatubo-scale peak (about 0.15); adjustable.
const SERIES_END_WARN_AOD = 0.005

# Schema version of the JLD2 file `load_aerosol_series` reads.
const AEROSOL_SERIES_FORMAT = 1

"""
    AerosolScenario(; eruptions=[], injections=[], series=nothing, tau_rise=0.21,
                    tau_decay=1.0, ssa=1.0, asymmetry=0.7, mu0=0.5)

Stratospheric aerosol add-on for `cfg.aerosol`. It composes with any
experiment and dims the shortwave per latitude. Optical depth follows two
linear reservoirs (rise timescale `tau_rise`, decay timescale `tau_decay`,
both in years; the defaults put the peak five months after injection and
decay with an e-folding time of one year, as for the tropical Pinatubo
aerosol); contributions from separate eruptions and injections add. The same
timescales apply to every injection class, although high-latitude eruptions
decay faster. `ssa`, `asymmetry` and `mu0` set the delta-Eddington
transmission.

`series` adds a published optical-depth record (see
[`load_aerosol_series`](@ref)), summed with the modelled eruptions. If the
record already contains an eruption you also list, it is counted twice. The
record is zero outside its first and last year, so a record that ends with a
non-zero value stops abruptly. Published records include the background
aerosol of quiet years (optical depth of about 0.002 to 0.009), which the
control run does not have, so a series run's anomaly contains that small
offset as well as the eruptions.

Applies to the scenario run only; the control run is never dimmed, so
anomalies are relative to an aerosol-free control. Not represented: the
longwave effect of the aerosol, its ozone chemistry, growth of particle size
over the eruption's lifetime (`ssa` and `asymmetry` are fixed), absorption in
the 8-12 micron window, and stratospheric heating.
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
end

function AerosolScenario(; eruptions=Eruption[], injections=SustainedInjection[],
        series::Union{Nothing,AerosolSeries}=nothing, tau_rise::Real=0.21,
        tau_decay::Real=1.0, ssa::Real=1.0, asymmetry::Real=0.7, mu0::Real=0.5)
    0 < tau_rise < tau_decay ||
        throw(ArgumentError("need 0 < tau_rise < tau_decay (years), got $tau_rise and $tau_decay"))
    0 <= ssa <= 1 || throw(ArgumentError("ssa must be in [0, 1], got $ssa"))
    0 <= asymmetry < 1 || throw(ArgumentError("asymmetry must be in [0, 1), got $asymmetry"))
    0 < mu0 <= 1 || throw(ArgumentError("mu0 must be in (0, 1], got $mu0"))
    return AerosolScenario(collect(Eruption, eruptions), collect(SustainedInjection, injections),
        series, Float64(tau_rise), Float64(tau_decay), Float64(ssa), Float64(asymmetry), Float64(mu0))
end

# - Two-reservoir chain (times in years) ──────────────────────────────────

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

Fills `aod` (length `ydim`) with the optical depth at decimal year `t`: the sum
over eruptions and sustained injections, each spread by its class profile.
Pure function of `t`; nothing persists between calls, and nothing is allocated.
"""
function aerosol_optical_depth!(aod::Vector{Float32}, sc::AerosolScenario, t::Real)
    fill!(aod, 0.0f0)
    tr, td = sc.tau_rise, sc.tau_decay
    inv_peak = 1 / chain_impulse(chain_peak_time(tr, td), tr, td)  # scales each pulse to a peak of 1
    for e in sc.eruptions
        a = e.peak_aod * inv_peak * chain_impulse(t - e.year, tr, td)
        a == 0.0 && continue
        aod .+= a .* _CLASS_PROFILES[e.cls]
    end
    for s in sc.injections
        u = chain_step(t - s.start_year, tr, td)
        isfinite(s.stop_year) && (u -= chain_step(t - s.stop_year, tr, td))
        u == 0.0 && continue
        aod .+= (s.target_aod * u) .* _CLASS_PROFILES[s.cls]
    end
    sc.series === nothing || _add_series!(aod, sc.series, t)
    return aod
end

# Adds the record at decimal year `t`, linear in time between records and zero
# outside the first and last year.
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

# - Delta-Eddington transmission ──────────────────────────────────────────
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
decimal year `t`, using `aod` (length `ydim`) as scratch for the optical depth
so nothing is allocated per timestep. Rows with zero optical depth are left
untouched, so the multiplier is bit-identical outside an eruption.
"""
function aerosol_transmission!(mult::Vector{Float32}, aod::Vector{Float32}, sc::AerosolScenario, t::Real)
    aerosol_optical_depth!(aod, sc, t)
    for j in 1:ydim
        tau = Float64(aod[j])
        tau == 0.0 && continue
        mult[j] *= Float32(delta_eddington(tau, sc.ssa, sc.asymmetry, sc.mu0)[2])
    end
    return mult
end
