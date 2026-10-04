# The dataset's file list, shared by the converter and the packager. It is the
# tools' own list, so a file can be added here before the model reads it.
# test/test_io.jl checks that every file the loaders read
# (GREBClimate.dataset_field_files) is in MODEL_FIELD_NAMES.

# Combined files that hold several fields. The converter writes them from
# special-case code, so they are not per-field entries below.
const FLUX_CORRECTION_NAMES = ("Tsurf_flux_correction", "vapour_flux_correction", "Tocean_flux_correction")
const COMBINED_FILE_NAMES = ("flux_corrections", "ipcc_scenarios", "solar_paleo",
                             "solar_eccentricity", "solar_obliquity",
                             "historical_emissions_population")

const MODEL_FIELD_NAMES = Set{String}([
    # static
    "global.topography", "greb.glaciers",
    # solar
    "solar_radiation.clim",
    # NCEP climatology (dataset=:ncep)
    "ncep.tsurf.1948-2007.clim", "ncep.zonal_wind.850hpa.clim",
    "ncep.meridional_wind.850hpa.clim", "ncep.atmospheric_humidity.clim",
    "ncep.soil_moisture.clim",
    # ERA-Interim climatology (dataset=:era; soil moisture falls back to NCEP)
    "erainterim.tsurf.1979-2015.clim", "erainterim.zonal_wind.850hpa.clim",
    "erainterim.meridional_wind.850hpa.clim", "erainterim.atmospheric_humidity.clim",
    # common to both datasets
    "isccp.cloud_cover.clim", "woce.ocean_mixed_layer_depth.clim", "Tocean.clim",
    "erainterim.omega.vertmean.clim", "erainterim.omega_std.vertmean.clim",
    "erainterim.windspeed.850hpa.clim",
    # CMIP5 RCP8.5 climate-change anomalies (load_cc_anomaly_jld2!)
    "cmip5.tsurf.rcp85.ensmean.forcing", "cmip5.zonal.wind.rcp85.ensmean.forcing",
    "cmip5.meridional.wind.rcp85.ensmean.forcing", "cmip5.windspeed.rcp85.ensmean.forcing",
    "cmip5.omega.rcp85.ensmean.forcing",
    # ENSO anomalies (load_enso_anomaly_jld2!), suffix elnino/lanina
    ("erainterim.$f.$s.forcing" for f in ("tsurf", "zonal.wind", "meridional.wind",
                                          "windspeed", "omega")
                                for s in ("elnino", "lanina"))...,
])
