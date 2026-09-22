"""Fetch an ERA5 monthly climatology and write GREB .bin files.

Writes flat 32-bit little-endian, no-header, Fortran-order (lon x lat x time)
binaries in the layout tools/convert_greb_to_jld2.jl reads. Field/level choices
were validated against the repo's existing erainterim.*.clim climatology.

Requires: cdsapi, netCDF4, numpy, and a configured ~/.cdsapirc
(https://cds.climate.copernicus.eu/how-to-api). Not a dependency of the
Julia package.

Usage:
    python -m pip install cdsapi netCDF4 numpy
    python tools/fetch_era5_data.py [--start-year 1991] [--end-year 2020] \\
        [--raw-dir era5_raw] [--out-dir era5_greb]
    julia --project=. tools/convert_greb_to_jld2.jl era5_greb greb_input_data_era5

Diagnostic fields (cloud_cover, t700, t1000, rh700) are fetched as raw netCDF
but not converted or written as GREB inputs; they support a separate
cloud-feedback experiment tracked in the vault, not this script.
"""

import argparse
import os

import cdsapi
import netCDF4 as nc
import numpy as np

GREB_YDIM = 48
GREB_NSTEP_YR = 730  # 365 days x 2 half-daily steps
GREB_LAT = np.array([-88.125 + 3.75 * i for i in range(GREB_YDIM)])  # cell centers, matches the repo's .ctl files

# name -> (CDS variable, netCDF variable name, dataset family, level or list of levels for a vertical mean, or None)
GREB_FIELDS = {
    "tsurf": ("2m_temperature", "t2m", "single-levels", None),
    "zonal_wind": ("u_component_of_wind", "u", "pressure-levels", 850),
    "meridional_wind": ("v_component_of_wind", "v", "pressure-levels", 850),
    "humidity": ("specific_humidity", "q", "pressure-levels", 1000),
    "omega": ("vertical_velocity", "w", "pressure-levels", [850, 700, 600, 500]),
}

# Fetched and converted for a separate cloud-feedback experiment; not written as GREB inputs.
DIAGNOSTIC_FIELDS = {
    "cloud_cover": ("total_cloud_cover", "tcc", "single-levels", None),
    "t700": ("temperature", "t", "pressure-levels", 700),
    "t1000": ("temperature", "t", "pressure-levels", 1000),
    "rh700": ("relative_humidity", "r", "pressure-levels", 700),
}

ALL_FIELDS = {**GREB_FIELDS, **DIAGNOSTIC_FIELDS}


def fetch_field(client, name, variable, family, level, start_year, end_year, raw_dir):
    """Download one field's monthly-means netCDF if not already present."""
    dataset = ("reanalysis-era5-single-levels-monthly-means" if family == "single-levels"
               else "reanalysis-era5-pressure-levels-monthly-means")
    request = {
        "product_type": "monthly_averaged_reanalysis",
        "variable": variable,
        "year": [str(y) for y in range(start_year, end_year + 1)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "time": "00:00",
        "grid": "3.75/3.75",
        "data_format": "netcdf",
    }
    if level is not None:
        levels = level if isinstance(level, list) else [level]
        request["pressure_level"] = [str(lv) for lv in levels]

    target = os.path.join(raw_dir, f"era5_{name}_raw.nc")
    if os.path.exists(target):
        print(f"{name}: {target} already exists, skipping fetch")
        return target
    print(f"Fetching {name} ({variable}, {start_year}-{end_year})...")
    client.retrieve(dataset, request, target)
    return target


def load_monthly_climatology(nc_path, ncvar, n_years):
    """Return monthly mean/std on the native grid, averaged over years and levels."""
    ds = nc.Dataset(nc_path)
    var = ds.variables[ncvar]
    dims = var.dimensions  # capture before close() - accessing it after close() raises
    data = np.asarray(var[:], dtype=np.float64)
    lat = np.asarray(ds.variables["latitude"][:])
    ds.close()

    if data.shape[0] != n_years * 12:
        raise ValueError(f"{nc_path}: expected {n_years * 12} time steps, got {data.shape[0]}")

    if "pressure_level" in dims:
        n_levels = data.shape[dims.index("pressure_level")]
        data = data.reshape(n_years, 12, n_levels, len(lat), -1)
        data = data.mean(axis=2)  # vertical mean across the requested levels
    else:
        data = data.reshape(n_years, 12, len(lat), -1)

    monthly_mean = data.mean(axis=0)          # (12, lat, lon)
    monthly_std = data.std(axis=0)            # (12, lat, lon) - interannual spread
    return monthly_mean, monthly_std, lat


def regrid_latitude(monthly, src_lat):
    """Interpolate latitude to GREB cell centers; longitude already matches."""
    order = np.argsort(src_lat)
    src_sorted = src_lat[order]
    n_month, _, n_lon = monthly.shape
    out = np.empty((n_month, GREB_YDIM, n_lon), dtype=np.float64)
    for m in range(n_month):
        for j in range(n_lon):
            out[m, :, j] = np.interp(GREB_LAT, src_sorted, monthly[m, order, j])
    return out


def monthly_to_annual_cycle(monthly):
    """Interpolate 12 monthly maps to 730 half-daily steps."""
    mid_month_doy = np.array([15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349], dtype=float)
    ext_doy = np.concatenate(([mid_month_doy[-1] - 365.0], mid_month_doy, [mid_month_doy[0] + 365.0]))
    _, n_lat, n_lon = monthly.shape
    daily = np.empty((365, n_lat, n_lon), dtype=np.float64)
    days = np.arange(1, 366, dtype=float)
    for i in range(n_lat):
        for j in range(n_lon):
            ext_vals = np.concatenate(([monthly[-1, i, j]], monthly[:, i, j], [monthly[0, i, j]]))
            daily[:, i, j] = np.interp(days, ext_doy, ext_vals)
    annual = np.empty((n_lat, n_lon, GREB_NSTEP_YR), dtype=np.float32)
    annual[:, :, 0::2] = daily.transpose(1, 2, 0)
    annual[:, :, 1::2] = daily.transpose(1, 2, 0)
    return annual  # (lat, lon, 730)


def to_greb_field(monthly, src_lat):
    """Convert monthly native-grid data to GREB (lon, lat, 730) float32."""
    regridded = regrid_latitude(monthly, src_lat)   # (12, GREB_YDIM, lon)
    annual = monthly_to_annual_cycle(regridded)      # (GREB_YDIM, lon, 730)
    return annual.transpose(1, 0, 2).astype(np.float32)  # (lon, lat, 730)


def write_bin(path, data):
    np.asarray(data, dtype=np.float32).flatten(order="F").tofile(path)
    print(f"Wrote {path}  shape={data.shape}  bytes={os.path.getsize(path)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--start-year", type=int, default=1991)
    ap.add_argument("--end-year", type=int, default=2020)
    ap.add_argument("--raw-dir", default="era5_raw", help="directory for downloaded netCDF files")
    ap.add_argument("--out-dir", default="era5_greb", help="directory for generated GREB .bin files")
    args = ap.parse_args()

    os.makedirs(args.raw_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    period = f"{args.start_year}-{args.end_year}"
    n_years = args.end_year - args.start_year + 1

    client = cdsapi.Client()
    raw_paths = {}
    for name, (variable, _ncvar, family, level) in ALL_FIELDS.items():
        raw_paths[name] = fetch_field(client, name, variable, family, level,
                                       args.start_year, args.end_year, args.raw_dir)

    fields_greb = {}
    omega_std_greb = None
    for name, (_variable, ncvar, _family, _level) in GREB_FIELDS.items():
        monthly_mean, monthly_std, lat = load_monthly_climatology(
            raw_paths[name], ncvar, n_years)
        fields_greb[name] = to_greb_field(monthly_mean, lat)
        if name == "omega":
            omega_std_greb = to_greb_field(monthly_std, lat)

    windspeed = np.sqrt(fields_greb["zonal_wind"] ** 2 + fields_greb["meridional_wind"] ** 2).astype(np.float32)

    write_bin(os.path.join(args.out_dir, f"era5.tsurf.{period}.clim.bin"), fields_greb["tsurf"])
    write_bin(os.path.join(args.out_dir, f"era5.zonal_wind.850hpa.{period}.clim.bin"), fields_greb["zonal_wind"])
    write_bin(os.path.join(args.out_dir, f"era5.meridional_wind.850hpa.{period}.clim.bin"), fields_greb["meridional_wind"])
    write_bin(os.path.join(args.out_dir, f"era5.atmospheric_humidity.{period}.clim.bin"), fields_greb["humidity"])
    write_bin(os.path.join(args.out_dir, f"era5.windspeed.850hpa.{period}.clim.bin"), windspeed)
    write_bin(os.path.join(args.out_dir, f"era5.omega.vertmean.{period}.clim.bin"), fields_greb["omega"])
    write_bin(os.path.join(args.out_dir, f"era5.omega_std.vertmean.{period}.clim.bin"), omega_std_greb)

    print("\nDone. Convert with:")
    print(f"  julia --project=. tools/convert_greb_to_jld2.jl {args.out_dir} <jld2-out-dir>")
    print("Wiring dataset=:era5 into GREBClimate.jl itself (flux-correction "
          "re-derivation, loading code) is a separate step - see the vault.")


if __name__ == "__main__":
    main()
