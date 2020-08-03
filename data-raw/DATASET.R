## code to prepare `DATASET` dataset goes here

# GW data
GWL <- readRDS("D:/PhD/HOAL/processed_data/piezometer/GWLabsLatest.RDS")
GWdepth <- readRDS("D:/PhD/HOAL/processed_data/piezometer/GWLdepthLatest.RDS")
GWtemp <- readRDS("D:/PhD/HOAL/processed_data/piezometer/GWLtempLatest.RDS")

usethis::use_data(GWL,
                  GWdepth,
                  GWtemp, overwrite = T)

# Discharge data
HOAL.Discharge <- readRDS("D:/PhD/HOAL/data_lovrenc/stream/ExportData/Q_all_2011-2017.RDS")

usethis::use_data(HOAL.Discharge, overwrite = T)

# Soil Moisture
HOAL.SM <- readRDS("D:/PhD/HOAL/data_lovrenc/soil moisture/SM_ISMN_13_18.RDS")
HOAL.SM.ASI <- readRDS("D:/PhD/HOAL/data_lovrenc/soil moisture/SM_ISMN_ASI_13_18.RDS")
HOAL.SM.Tot <- readRDS("D:/PhD/HOAL/data_lovrenc/soil moisture/SM_ISMN_tot_13_18.RDS")

usethis::use_data(HOAL.SM,
                  HOAL.SM.ASI,
                  HOAL.SM.Tot, overwrite = T)

# Weather station data
HOAL.Air.Press <- readRDS("D:/PhD/HOAL/processed_data/weather_station/Atmos_pressure_12-18.RDS")
HOAL.Air.Temp <- readRDS("D:/PhD/HOAL/processed_data/weather_station/T_RH_dat_quality_12-18.RDS")
HOAL.Soil.Temp <- readRDS("D:/PhD/HOAL/processed_data/weather_station/Soil_temp_12-18.RDS")
HOAL.Rad <- readRDS("D:/PhD/HOAL/processed_data/weather_station/Rad_data_quality_12-18.RDS")

usethis::use_data(HOAL.Air.Press,
                  HOAL.Air.Temp,
                  HOAL.Soil.Temp,
                  HOAL.Rad, overwrite = T)

# Rainfall
HOAL.Rain <- readRDS("D:/PhD/HOAL/data_lovrenc/rainfall/ExportData/HOAL_rain_10-18.RDS")
HOAL.Rain.I <- readRDS("D:/PhD/HOAL/data_lovrenc/rainfall/ExportData/HOAL_intensity_10-18.RDS")
HOAL.Rainfall <- readRDS("D:/PhD/HOAL/data_lovrenc/rainfall/ExportData/HOAL_rainfall_10-18.RDS")
HOAL.Rain.Mean5 <- readRDS("D:/PhD/HOAL/data_lovrenc/rainfall/ExportData/HOAL_rainfall_10-18_mean_5min.RDS")

usethis::use_data(HOAL.Rain,
                  HOAL.Rainfall,
                  HOAL.Rain.I,
                  HOAL.Rain.Mean5, overwrite = T)

## GIS data
HOAL.CRS <- sp::CRS("+init=epsg:31256")
HOAL.boundary <- rgdal::readOGR("D:/PhD/HOAL/GIS/catchment boundary",
                      "catchment_boundary")
HOAL.GW.stations <- rgdal::readOGR("D:/PhD/HOAL/GIS/piezometers", "GW_stations")
HOAL.Q.stations <- rgdal::readOGR("D:/PhD/HOAL/GIS/discharge", "discharge_monitoring")
HOAL.stream <- rgdal::readOGR("D:/PhD/HOAL/GIS/stream", "stream_dhm1")
HOAL.contour <- rgdal::readOGR("D:/PhD/HOAL/data_lovrenc/GIS/DEM","contours_05m")
HOAL.ortophoto <- raster::stack("D:/PhD/HOAL/GIS/orthofoto/Orthofoto_Petzenkirchen_GK_East_31256.tif")

usethis::use_data(HOAL.CRS,
                  HOAL.boundary,
                  HOAL.GW.stations,
                  HOAL.Q.stations,
                  HOAL.stream,
                  HOAL.contour,
                  HOAL.ortophoto, overwrite = T)
