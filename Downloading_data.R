####################################
######### DOWNLOADING DATA #########
####################################

# Downloading data from three reanalysis databases:
# NCEP/NCAR Reanalysis 1, NCEP/DOE Reanalysis 2 and ECMWF ERA5.

# IMPORTANT: To download data from ERA5, registration is required at the Copernicus Climate Data Store:
# https://cds.climate.copernicus.eu
# UserID and API key shall be given in the wf_set_key function in this script.


####################################################
######### DOWNLOADING NCEP REANALYSIS DATA #########
####################################################

######### Input parameters #########
# year            # time interval by choosing first and last years, e.g. c(1981,2020), default is 2020
# month           # list of months, e.g. c(1,3), 1:2, 1, default is 3 (March)
# latitude        # borders of the requested area: latitudes, e.g. c(0,90), default is 47.5dg
# longitude       # borders of the requested area: longitudes, e.g. c(-177.5,180), default is 20dg
# rean            # choosing database: FALSE - NCEP/NCAR Reanalysis 1, TRUE - NCEP/DOE Reanalysis 2
                  # default is FALSE
# variable        # choosing variable, e.g air temperature is "air", default is geopotential: "hgt"
# variable_long   # name of the variable in the file title, e.g. "geopot500"
# height          # pressure level given in hPa, e.g. 1000, default is 500

NCEPdownload <- function(year=c(2020,2020), month=3,
                         latitude=c(47.5,47.5), longitude=c(20,20), rean=FALSE, 
                         variable="hgt", variable_long="geopot500", height=500){

  # Data will be stored in the object vars:
  vars <- c(NA)
  
  ######### Download data #########
  vars <- NCEP.gather(variable=variable,
                      level=height,
                      months.minmax=month,
                      years.minmax=year,
                      lat.southnorth=latitude,
                      lon.westeast= longitude,
                      reanalysis2=rean,
                      return.units=TRUE,
                      status.bar=TRUE)
    saveRDS(object=vars, file=paste0(variable_long,"_NCEP","_",year[1],"-",year[2],".rds"))
}


####################################################
######### DOWNLOADING ERA5 REANALYSIS DATA #########
####################################################

######### Input parameters #########
# Identifier at Copernicus:
# user           # userID
# API            # APIkey
# year           # time interval by choosing first and last years, e.g. c(1981,2010), default is c(2020,2020)
# month          # list of months, e.g. c("01","03"), default is "03" (March)
# day            # list of days, e.g. as.character(1:31), default is "01"
# hour           # list of hours, e.g. c("00:00","06:00","12:00","18:00"), default is "12:00"
# lon1           # border of the requested area: west, default is 18.8
# lon2           # border of the requested area: east, default is 18.8
# lat1           # border of the requested area: south, default is 47.3
# lat2           # border of the requested area: north, default is 47.3
# variable       # choosing variable, e.g. "2m_temperature", default is "geopotential"
# variable_long  # name of the variable in the file title, default is "geopot500"
# rean           # choosing database: "reanalysis-era5-pressure-levels",
                 # or "reanalysis-era5-single-levels" to get surface level data
# resolution     # resolution of the grid in degree, e.g. "2.5/2.5", default is "0.1/0.1"
# height         # pressure level given in hPa, e.g. 500, set it as 0 to get surface data
# timeOut        # how long to wait to download data in seconds, e.g. 1800
# fileDir        # folder in which data will be downloaded, default is the working directory

ERA5download <- function(year=c(2020,2020), month=03, day=01, hour="12:00", lat1=47.3, lat2=47.3,
                         lon1=18.8, lon2=18.8, rean="reanalysis-era5-pressure-levels",
                         variable="geopotential", variable_long="geopot500", resolution="0.1/0.1",
                         height=500, timeOut=1800, fileDir=getwd()){
 # Identifying user:
  userID <- wf_set_key(user = user,
                       key = API,
                       service = "cds")
  
  ######### Download data #########
    if(height > 0 & rean == "reanalysis-era5-pressure-levels") { # downloading to pressure level
      request <- list(
        "dataset_short_name" = rean,
        "class" = "ea",
        "expver" = "1",
        "stream" = "oper",
        "product_type" = "reanalysis",
        "pressure_level" = height,
        "variable" = variable,
        "year" = as.character(year[1]:year[2]),
        "month" = as.character(month),
        "day" = as.character(day),
        "time" = hour,
        "area" = paste0(lat1,"/",lon1,"/",lat2,"/",lon2), # N, W, S, E
        "grid" = resolution,
        "format" = "netcdf",
        "target" = paste0(variable_long,"_ERA5","_",year[1],"-",year[2],".nc"))
      
    file <- wf_request(user     = userID,
                       request  = request,  
                       transfer = TRUE,
                       path     = fileDir,
                       time_out = timeOut)
    }
    if(height == 0 & rean == "reanalysis-era5-single-levels") {# downloading to surface level
      request <- list(
        "dataset_short_name" = rean,
        "class" = "ea",
        "expver" = "1",
        "stream" = "oper",
        "product_type" = "reanalysis",
        "variable" = variable,
        "year" = as.character(year[1]:year[2]),
        "month" = month,
        "day" = day,
        "time" = hour,
        "area" = paste0(lat1,"/",lon1,"/",lat2,"/",lon2), # N, W, S, E
        "grid" = resolution,
        "format" = "netcdf",
        "target" = paste0(variable_long,"_ERA5","_",year[1],"-",year[2],".nc"))
      
      file <- wf_request(user     = userID,
                         request  = request,  
                         transfer = TRUE,
                         path     = fileDir,
                         time_out = timeOut)
    
  }
}
