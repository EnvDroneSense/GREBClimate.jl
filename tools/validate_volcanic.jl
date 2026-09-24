### Volcanic validation: model response to Agung, El Chichon and Pinatubo against GISTEMP ###
#
# MAINTAINER TOOL - not part of the package. Needs the local dataset, the
# converted aerosol records and the validation data in Data/:
#
#   Data/aerosol/aerosol_Sato-Lacis.jld2, aerosol_GloSSAC.v2.24.jld2
#       tools/convert_aerosol_to_jld2.jl on the GISS tau_reff_*.nc files
#   Data/validation/GLB.Ts+dSST.csv
#       https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.csv
#   Data/validation/mei_v1_table.html
#       https://psl.noaa.gov/enso/mei.old/table.html (MEI, 1950-2018)
#
# 1. Runs :full_model with and without each aerosol input; the difference in
#    global-mean surface temperature is the model's volcanic response.
# 2. Regresses the observed monthly anomaly on a trend, the lagged MEI, a
#    two-harmonic annual cycle and the model response in a 5-year window after
#    each eruption (Foster and Rahmstorf 2011, with the model response in
#    place of optical depth). A coefficient near 1 means the model has the
#    right magnitude for that event.
# 3. Writes the time series (observations with trend, ENSO and annual cycle
#    removed, and the model response) to Data/validation/volcanic_timeseries.csv.
#
# Usage:
#   julia --project=. tools/validate_volcanic.jl

using GREBClimate
using LinearAlgebra
using Printf

const REPO = normpath(joinpath(@__DIR__, ".."))
const DATA_DIR = something(greb_data_dir(; allow_download=false),
                           joinpath(REPO, "greb_input_data"))
const AEROSOL_DIR = joinpath(REPO, "Data", "aerosol")
const VALIDATION_DIR = joinpath(REPO, "Data", "validation")

const FIRST_YEAR = 1950                  # the scenario run starts here
const LAST_YEAR = 2018                   # end of MEI v1
const RUN = RunSpec(flux=3, ctrl=5, scnr=LAST_YEAR - FIRST_YEAR + 1)
const WINDOW = 60                        # months attributed to each eruption
const BASELINE = 12                      # months before the eruption used as its baseline

# (name, year, month) of each eruption
const ERUPTIONS = (("Agung", 1963, 3), ("El Chichon", 1982, 4), ("Pinatubo", 1991, 6))

month_index(year, month) = 12 * (year - FIRST_YEAR) + month
decimal_year(m) = FIRST_YEAR + (m - 0.5) / 12

# ── Model runs ───────────────────────────────────────────────────────────

function global_mean_ts(result)
    w = cosd.(Float64.(GREBClimate.lat_grid))
    return [sum(rec.Ts .* reshape(w, 1, :)) / (xdim * sum(w)) for rec in result.scnr]
end

function run_model(fields, aerosol)
    cfg = create_experiment_config(:full_model)
    cfg.aerosol = aerosol
    result = redirect_stdout(devnull) do
        greb_model!(RUN, cfg; jld2_dir=DATA_DIR, fields=deepcopy(fields))
    end
    return global_mean_ts(result)
end

# ── Observations ─────────────────────────────────────────────────────────

# GISTEMP monthly global means (deg C, 1951-1980 base), one per model month.
function load_gistemp(path)
    y = fill(NaN, 12 * (LAST_YEAR - FIRST_YEAR + 1))
    for line in eachline(path)
        cols = split(line, ',')
        yr = tryparse(Int, cols[1])
        (yr === nothing || !(FIRST_YEAR <= yr <= LAST_YEAR)) && continue
        for mo in 1:12
            v = tryparse(Float64, cols[mo + 1])
            v === nothing || (y[month_index(yr, mo)] = v)
        end
    end
    return y
end

# MEI v1 bimonthly values (DECJAN .. NOVDEC), each assigned to the second month.
function load_mei(path)
    e = fill(NaN, 12 * (LAST_YEAR - FIRST_YEAR + 1))
    for line in eachline(path)
        cols = split(strip(line))
        isempty(cols) && continue
        yr = tryparse(Int, cols[1])
        (yr === nothing || !(FIRST_YEAR <= yr <= LAST_YEAR)) && continue
        for (mo, s) in enumerate(cols[2:min(end, 13)])
            v = tryparse(Float64, s)
            v === nothing || (e[month_index(yr, mo)] = v)
        end
    end
    return e
end

# ── Regression ───────────────────────────────────────────────────────────

# The model response in the window after each eruption, relative to the mean
# of the months before it (removes the background aerosol), zero elsewhere.
function event_regressors(response, events)
    n = length(response)
    M = zeros(n, length(events))
    for (k, (_, yr, mo)) in enumerate(events)
        m0 = month_index(yr, mo)
        base = sum(response[m0 - BASELINE:m0 - 1]) / BASELINE
        for m in m0:min(n, m0 + WINDOW - 1)
            M[m, k] = response[m] - base
        end
    end
    return M
end

# OLS fit of y over months `rows` with MEI lag `lag`. Returns the coefficients,
# their standard errors inflated for AR(1) residuals, the residual sum of
# squares, the design matrix and the number of event columns.
function fit(y, mei, M, rows, lag; quadratic=false)
    t = [decimal_year(m) - 1980 for m in rows]
    phase = [2pi * (m - 1) / 12 for m in rows]
    X = hcat(ones(length(rows)), t, (quadratic ? [t .^ 2] : [])...,
             mei[rows .- lag], cos.(phase), sin.(phase), cos.(2 .* phase), sin.(2 .* phase),
             M[rows, :])
    yy = y[rows]
    beta = X \ yy
    r = yy - X * beta
    rho = dot(r[1:end-1], r[2:end]) / dot(r, r)
    sigma2 = dot(r, r) / (length(rows) - size(X, 2))
    se = sqrt.(sigma2 .* diag(inv(X' * X)) .* (1 + rho) / (1 - rho))
    return (beta=beta, se=se, rss=dot(r, r), X=X, rho=rho)
end

function report(label, y, mei, M, rows, events)
    ne = length(events)
    rss = [fit(y, mei, M, rows, lag).rss for lag in 0:24]
    best = argmin(rss) - 1
    f = fit(y, mei, M, rows, best)
    fq = fit(y, mei, M, rows, best; quadratic=true)
    ranges = [extrema(fit(y, mei, M, rows, lag).beta[end - ne + k] for lag in 0:6) for k in 1:ne]
    @printf("\n%s: %d-%d, best MEI lag %d months, residual AR(1) %.2f\n",
            label, floor(Int, decimal_year(first(rows))), floor(Int, decimal_year(last(rows))),
            best, f.rho)
    println("  event         beta   +/- 2 se    lags 0-6        quadratic trend")
    for k in 1:ne
        i = length(f.beta) - ne + k
        @printf("  %-12s %5.2f   %5.2f      %5.2f to %5.2f   %5.2f\n",
                events[k][1], f.beta[i], 2 * f.se[i], ranges[k]..., fq.beta[end - ne + k])
    end
    return (lag=best, fit=f)
end

# ── Main ─────────────────────────────────────────────────────────────────

function main()
    fields = load_greb_jld2!(DATA_DIR; dataset=:ncep)
    sato = load_aerosol_series(joinpath(AEROSOL_DIR, "aerosol_Sato-Lacis.jld2"))
    glossac = load_aerosol_series(joinpath(AEROSOL_DIR, "aerosol_GloSSAC.v2.24.jld2"))

    println("running the model (", RUN.scnr, "-year scenario, 6 runs)...")
    base = run_model(fields, nothing)
    resp = Dict(
        "sato" => run_model(fields, AerosolScenario(series=sato)) - base,
        "glossac" => run_model(fields, AerosolScenario(series=glossac)) - base)
    for tg in (7, 9, 11)
        resp["mass$tg"] = run_model(fields,
            AerosolScenario(eruptions=[Eruption(1991.45, :tropical; tg_s=tg)])) - base
    end

    # Mass mode: Pinatubo from the injected sulfur, 7 to 11 Tg S (Toohey and Sigl 2017).
    println("\nPinatubo, mass mode (3-month running mean of the global-mean response)")
    m0 = month_index(1991, 6)
    for tg in (7, 9, 11)
        r = resp["mass$tg"]
        s = [sum(r[m-1:m+1]) / 3 for m in m0:m0 + WINDOW]
        i = argmin(s)
        rec = findfirst(x -> x > s[i] / exp(1), s[i:end])
        @printf("  %2d Tg S: peak %.2f K, %d months after the eruption; 1/e recovery %s\n",
                tg, s[i], i - 1, rec === nothing ? "beyond 5 years" : "after $(rec - 1) more months")
    end

    y = load_gistemp(joinpath(VALIDATION_DIR, "GLB.Ts+dSST.csv"))
    mei = load_mei(joinpath(VALIDATION_DIR, "mei_v1_table.html"))
    last_m = findlast(!isnan, mei)                       # MEI v1 stops in 2018

    M_sato = event_regressors(resp["sato"], ERUPTIONS)
    rows_sato = month_index(1952, 1):month_index(2012, 12)   # lag room; Sato ends 2012
    s = report("Sato-Lacis", y, mei, M_sato, rows_sato, ERUPTIONS)

    later = ERUPTIONS[2:3]                               # GloSSAC starts in 1979
    M_glossac = event_regressors(resp["glossac"], later)
    report("GloSSAC v2.24", y, mei, M_glossac, month_index(1981, 1):last_m, later)

    # Time series: observations with everything but the eruption terms removed.
    f = s.fit
    ne = length(ERUPTIONS)
    other = f.X[:, 1:end - ne] * f.beta[1:end - ne]
    path = joinpath(VALIDATION_DIR, "volcanic_timeseries.csv")
    open(path, "w") do io
        println(io, "year,observed_adjusted,model_sato,model_glossac,model_sato_fitted")
        for (i, m) in enumerate(rows_sato)
            fitted = dot(M_sato[m, :], f.beta[end - ne + 1:end])
            @printf(io, "%.4f,%.4f,%.4f,%.4f,%.4f\n", decimal_year(m), y[m] - other[i],
                    resp["sato"][m], resp["glossac"][m], fitted)
        end
    end
    println("\nwrote ", path)
end

main()
