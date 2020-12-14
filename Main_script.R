#########################################################
###################### DESCRIPTION ######################
#########################################################

### A) This script downloads gridded time series of geopotential data from three reanalyses
###    to detect atmospheric teleconnections by using the functions in Downloading_data.R.

### B) Monthly geopotential height anoamlies are calculated on which correlation analysis 
###    (i.e. teleconnectivity method) (1) and principal component analysis (PCA) (2) are executed
###    by using the function in Data_procession_and_computation.R.

### C) Based on (1) strongest negative Pearson correlations are computed in each grid cell
###    and potential action centers (PotACs) are determined. (PotACs are grid cells that
###    are associated with correlations which are in one-to-one correspondence with each other.)
###    The first principal component (PC1) time series is examined.

### D) The fields of strongest negative correlations and the first empirical orthogonal
###    function (EOF1) are visualized on maps. Customized colorbar is made by applying
###    the function in "Creating_colorbar.R". Barplots of the PC1 times series are created.

### IMPORTANT: Data of NCEP-NCAR Reanalysis 1, NCEP-DOE Reanalysis 2 and ERA5 are used.
### To download data from ERA5, registration is required at the Copernicus Climate Data Store:
### https://cds.climate.copernicus.eu
### UserID and API key shall be given in the wf_set_key function in this script.
### The scripts Downloading_data.R, Data_procession_and_computation.R and Creating_colorbar.R
### shall be stored in the working directory or their path shall be given below.

### Data are downloaded in netCDF file in case of ERA5, and in rds format in case of NCEP.
### Results of correlation analysis and principal component analysis are stored in rds files.
### Images are created in png format.
###

### Data will be downloaded and processed in the working directory which can be set here:
setwd("...") # e.g. setwd("C:/Stat_analysis")

### CITATIONS OF DATABASES ###

# Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Munoz-Sabater, J., Nicolas, J., 
# Peubey, C., Radu, R., Schepers, D., Simmons, A., et al., 2020: The ERA5 global reanalysis. Quarterly
# Journal of the Royal Meteorological Society, 146, 1999–2049. https://doi.org/10.1002/qj.3803

# Kalnay, E., Kanamitsu, M., Kistler, R., Collins, W., Deaven, D., Gandin, L., Iredell, M., Saha, S.,
# White, G., Woollen, J., et al., 1996: The NCEP/NCAR 40-Year Reanalysis Project. Bulletin of the
# American Meteorological Society, 77(3), 437–472.
# https://doi.org/10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2

# Kanamitsu, M., Ebisuzaki, W., Woollen, J., Yang S-K., Hnilo, J.J., Fiorino, M., Potter, G.L., 2002:
# NCEP-DOE AMIP-II Reanalysis (R-2). Bulletin of the American Meteorological Society, 83(11), 1631–1644.
# https://doi.org/10.1175/BAMS-83-11-1631


############################################################
############## REQUIRED PACKAGES AND SCRIPTS ###############
############################################################

### Load required packages:
library(RNCEP)          # NCEP.gather: to download NCEP reanalysis data
library(ecmwfr)         # wf_set_key, wf_request: to download ERA5 reanalysis data
library(ncdf4)          # nc_open, ncvar_get, nc_close: to handle netCDF file
library(maps)           # map: to fit country borders to the maps
library(RColorBrewer)   # colorRampPalette, brewer.pal: to create colorscale
library(Hmisc)          # rcorr: to compute correlations between the gridded time series
library(magick)         # image_read, image_trim: to crop margins off

### Load required scripts:
source(paste("Downloading_data.R", sep="/"))                # NCEPdownload, ERA5download
source(paste("Data_procession_and_computation.R", sep="/")) # corr_pca_monthly_anom
source(paste("Creating_colorbar.R", sep="/"))               # colorbar_horizontal


#############################################################
#################### A) DOWNLOADING DATA ####################
#############################################################

### Downloading six hourly 500 hPa geopotential height data from three reanalysis databases.

### Input parameters of NCEPdownload function to get NCEP/NCAR Reanalysis 1 and NCEP/DOE Reanalysis 2 data:
year <- c(1981,2020)            # choosing time interval in years
month <- 1                      # giving the list of months
latitude <- c(0,90)             # borders of the requested area: latitudes
longitude <- c(-177.5,180)      # borders of the requested area: longitudes
rean <- FALSE                   # choosing database: FALSE - NCEP/NCAR Reanalysis 1, TRUE - NCEP/DOE Reanalysis 2
variable <- "hgt"               # choosing variable: geopotential height
variable_long <- "geopot500"    # name of the variable in the file title
height <- 500                   # pressure level given in hPa

# Download data:
NCEPdownload(year=year, month=month, latitude=latitude, longitude=longitude,
             rean=FALSE, variable=variable, variable_long=variable_long, height=height)

### Input parameter for ERA5download function to get ERA5 data:
# Identifier at Copernicus:
user <- ...
API  <- ...

year <- c(1981,2010)                          # choosing time interval in years
month <- "01"                                 # giving the list of months
day <- as.character(1:31)                     # giving the list of days
hour <- c("00:00","06:00","12:00","18:00")    # giving the list of hours
lon1 <- -177.5                                # border of the requested area: west
lon2 <- 180                                   # border of the requested area: east
lat1 <- 0                                     # border of the requested area: south
lat2 <- 90                                    # border of the requested area: north
rean <- "reanalysis-era5-pressure-levels"     # choosing the ERA5 reanalysis database
variable <- "geopotential"                    # choosing variable: geopotential height
variable_long <- "geopot500"                  # name of the variable in the file title
resolution <- "2.5/2.5"                       # resolution of the grid in degree
height <- 500                                 # pressure level given in hPa
timeOut <- 3600                               # how long to wait to download data
output_path <- getwd()                        # choosing downloading folder

# Download data:
ERA5download(year=year, month=month, day=day, hour=hour, lat1=lat1, lat2=lat2,
             lon1=lon1, lon2=lon2, rean=rean, variable=variable, variable_long=variable_long,
             resolution=resolution, height=height, timeOut=timeOut, fileDir=output_path)


#############################################################
############ B) DATA PREPROCESSION & COMPUTATION ############
#############################################################

### Input parameters:
input_file <- "geopot500_ERA5_1981-2020.nc" # "geopot500_ERA5_1981-2020.nc" # "geopot500_NCEP1_1981-2020.rds" # "geopot500_NCEP2_1981-2020.rds"

# Area, resolution (in degree) and time interval of the examination: North-Atlantic region
lon1A     <- -75
lon2A     <- 40
lat1A     <- 20
lat2A     <- 80
res       <- 2.5
startDate <- "1981-01-01 00:00:00"
endDate   <- "2010-01-31 00:00:00"

# Names of datasets:
id <- c("ERA5_81","ERA5_86","ERA5_91",  "NCEP1_81","NCEP1_86","NCEP1_91",  "NCEP2_81","NCEP2_86","NCEP2_91")

# Number of Monte Carlo simulation: to carry out permutation test on the correlation field:
experience_nr <- 1000

# Examined database:
# "NCEP" if NCEP-NCAR Reanalysis 1 or NCEP-DOE Reanalysis 2 is examined. "ERA5" if ERA5 is examined.
dataset <- "ERA5"

# Obtaining strongest negative correlation fields and do PCA from monthly averages:
corr_pca_monthly_anom(input_path=getwd(), input_file=input_file, dataset=dataset,
                     lon1=lon1A, lon2=lon2A, lat1=lat1A, lat2=lat2A, experience_nr=exp_nr,
                     start_date=startDate, end_date=endDate, output_path=getwd(),
                     output_file="absmincor_NCEP2_1981-2010_Jan.rds",
                     output_file2="pca_NCEP2_1981-2010_Jan.rds")

# Set interval in which correlations are examined:
cor_sign_min <- -0.65 # Based on permutation test, carried out by corr_pca_monthly_anom
                      # (distribution of randomly generated correlations are printed on the screen)
cor_sign_max <- -0.85 
step <- 0.025

# Loading the strongest negative correlation fields and the results of PCA:
files_abs <- Sys.glob(file.path("absmincor*Jan.rds"))
files_pca <- Sys.glob(file.path("pca*Jan.rds"))

res_list_abs <- list(NA)
for (mod in 1:length(files_abs)) {
  res_list_abs[[mod]] <- readRDS(files_abs[mod])
  res_list_abs[[mod]][res_list_abs[[mod]] > cor_sign_min] <- NA
}

res_list_pca <- list(NA)
for (mod in 1:length(files_pca)) {
  res_list_pca[[mod]] <- readRDS(files_pca[mod])
}


##########################################
############ C) DATA ANALYSIS ############
##########################################

################## FINDING POTENTIAL ACTION CENTERS (PotACs) ##################

coord_mod_list <- list(NA) # Empty list in which coordinates of PotACs will be stored.
lon <- seq(lon1A,lon2A,res) # Sequence of longitudes
lat <- seq(lat1A,lat2A,res) # Sequence of longitudes

potAC_with_stongest_corr <- data.frame(NA)
for (mod in 1:length(files_abs)) {
  
  coord_list <- list(NA)
  absmincor_vector <- na.omit(as.numeric(res_list_abs[[mod]]))
  last <- length(absmincor_vector[duplicated(absmincor_vector)])
  coord_matr <- array(NA, dim=c(last,5)) # Correlations which appear twice (i.e. PotACs!)
  
  for (vector in 1:last) {
    coord_list[[vector]] <- which(res_list_abs[[mod]] == absmincor_vector[duplicated(absmincor_vector)][vector], TRUE)

    coord_matr[vector,] <- cbind(
      lon[coord_list[[vector]][1,1]], lat[coord_list[[vector]][1,2]],
      lon[coord_list[[vector]][2,1]], lat[coord_list[[vector]][2,2]], 
      absmincor_vector[duplicated(absmincor_vector)][vector]
    )
  }
  
  try(if(nrow(coord_matr)!=last) stop("At least one PotAC is lost!!"))
  coord_matr <- as.data.frame(coord_matr)
  colnames(coord_matr) <- c("lon1","lat1","lon2","lat2","corr")

  coord_matr <- coord_matr[order(coord_matr$corr, decreasing=FALSE),]
  coord_mod_list[[mod]] <- coord_matr
  
  # Determining the strongest correlation and its coordinated based on each dataset:
  for (k in 1:ncol(coord_matr)) (potAC_with_stongest_corr[mod,k] <- coord_mod_list[[mod]][1,k])
}

colnames(potAC_with_stongest_corr) <- c("lon1","lat1","lon2","lat2","corr")
potAC_with_stongest_corr <- cbind(id, potAC_with_stongest_corr)


################## COMPARISON OF PC1 TIME SERIES ##################

### COMPARISON OF THE PC1 TIME SERIES ###

# Extract PC1 time series:
pc1_time_series <- list(NA)
for (mod in 1:length(files_pca)) ( pc1_time_series[[mod]] <- res_list_pca[[mod]]$x[,1] )

try(if(length(unlist(pc1_time_series)) %% length(files_pca) != 0) 
  stop("Time series do not have the same length!"))

# Order PC1 time series into a matrix.
#y Rows are time steps, columns are ERA5 (1981-2010), ERA5 (1986-2015), ..., NCEP-DOE (1991-2020)
pc1 <- matrix(unlist(pc1_time_series), nrow=length(unlist(pc1_time_series))/length(files_pca),
                                       ncol=length(files_pca))

if(!exists(id)) {
  colnames(pc1) <- id
  rownames(pc1) <- id
}

# Correlation matrix:
round(rcorr(pc1)$r, digits=3)

# If the 1st, 4th, 7th columns contain times series of the first examined period,
#    the 2nd, 5th, 8th columns contain times series of the first examined period and
#    the 3rd, 6th, 9th columns contain times series of the first examined period:
round(rcorr(pc1[,c(1,4,7)])$r, digits=3)
round(rcorr(pc1[,c(2,5,8)])$r, digits=3)
round(rcorr(pc1[,c(3,6,9)])$r, digits=3)


##########################################
############ D) VISUALIZATION ############
##########################################

######### Maps of strongest negative correlation fields and EOF1 #########

# Choosing the examined field:
examined_field <- "eof1" # "eof1" or "absmincor"

# Axis, labels and colorbars, plot name:
x <- c("","-60°","-40°","-20°","0°","20°","40°")
y <- c("20°","40°","60°","80°")
periods <- paste0("                         1981-2010                         ", 
                  "                   1986-2015                         ",
                  "                   1991-2020                         ")
reans <- c("ERA5","","","NCEP/NCAR R1","","","NCEP/DOE R2","","")
reans2 <- c("ERA5","NCEP/NCAR R1","NCEP/DOE R2")

# Range to create colorbars:
print(range(unlist(res_list_abs), na.rm=TRUE))
eof1_list <- list(NA)
for(mod in 1:length(res_list_pca)) (eof1_list[[mod]] <- res_list_pca[[mod]]$rotation[,1])
range(unlist(eof1_list))

# Axis, labels, colorbars, plot name:
if(examined_field == "absmincor") {
  plotname <- "absmincors.png"
  brks <- seq(cor_sign_max,cor_sign_min,step) # only significant corrs
  colors <- rev(colorRampPalette(brewer.pal(9,"Blues"))(length(brks)-1))
  colors_AC <- rev(colorRampPalette(brewer.pal(9,"Reds"))(length(brks)+2)) # Three more colors is added than requested.
  colors_AC <- colors_AC[1:(length(colors_AC)-3)] # The three brightest shades are omitted.
  res_list <- res_list_abs
}

if(examined_field == "eof1") {  
  plotname <- "eof1.png"
  plotname2 <- "pc1.png"
  brks <- seq(-0.06,0.06,0.01)
  colors <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(brks)-1))
  for (mod in 1:length(res_list_pca)) {
    res_list[[mod]] <- matrix(res_list_pca[[mod]]$rotation[,1], nrow=length(lon), ncol=length(lat))
  }
}
  
# Plotting data: strongest negative correlations fields with PotACs or EOF1 with explained variances
png(plotname, units="in", width=30, height=24, res=300, pointsize=54)
  par(omi=c(0,0,0,0)) # configure outer margin
  par(mfrow=c(4,3),
      mar=c(1.2,2.7,1.5,0.2)) # configure margin of each plot
  plot(x=1, y=1, pch="", ann=FALSE, axes=FALSE)
  par(new=TRUE)
  title(periods, outer=TRUE, line=-1 ,adj=0)
  par(new=TRUE)
  for (i in 1:length(res_list)) {
    par(mgp = c(0,0.3,0)) # tick labels are moved closer to the axes
    image(x=lon, y=lat, z=res_list[[i]], col=colors,
          breaks=brks, ann=FALSE, xaxt="n", yaxt="n")
    title(ylab=reans[i], cex.lab=1.3, font=1, line=1.7)
    abline(h=seq(-80,80,10), v=seq(-180,160,10), lty=2)
    map("world", interior=FALSE, add=TRUE)
    map("world", regions="Hungary", interior=FALSE, add=TRUE)
    
    # Associating colors to PotACs in case of strongest negative correlations:
    if (examined_field == "absmincor") {
      cols_of_AC <- c(NA)
      for (k in 1:nrow(coord_mod_list[[i]])){
        for (j in 1:(length(brks)-1)) {
          if(coord_mod_list[[i]]$cor[k] >= brks[j] & coord_mod_list[[i]]$cor[k]<brks[j+1]) {
            cols_of_AC[k] <- colors_AC[j]
           }
        }
      }
    }
      
    points(x=c(coord_mod_list[[i]]$lon1, coord_mod_list[[i]]$lon2),
           y=c(coord_mod_list[[i]]$lat1, coord_mod_list[[i]]$lat2), pch=16, col=cols_of_AC)
    arrows(x0=coord_mod_list[[i]]$lon1, y0=coord_mod_list[[i]]$lat1,
           x1=coord_mod_list[[i]]$lon2, y1=coord_mod_list[[i]]$lat2,
           lwd=2, col=cols_of_AC, code=0)
    
    if (examined_field == "eof1") {
      mtext(text=paste0(round(res_list_pca[[i]]$sdev[1], digits=0), " %"), side=3, adj=1, cex=0.7)
    }
    axis(1, at=seq(-80,40,20), labels=x)
    axis(1, at=seq(-80,40,10), labels=FALSE, tck=-0.02)
    axis(2, at=seq(20,80,20), labels=y, las=2)
    axis(2, at=seq(20,80,10), labels=FALSE, tck=-0.02)
    axis(3, at=c(-180,180), tick=TRUE, labels=FALSE, tck=0.01)
    axis(c(2,4), at=c(-90,90), tick=TRUE, labels=FALSE, tck=0.01)
  }
  if (examined_field == "absmincor") (par(fig=c(0.01,0.5,0.15,0.28), new=TRUE))
  if (examined_field == "eof1") (par(fig=c(0.12,0.88,0.16,0.31), new=TRUE))
  colorbar_horizontal(col=colors, brks=brks, height=4) # colorbar of absmincors in the 4th row
  if (examined_field == "absmincor") (par(fig=c(0.5,0.99,0.15,0.28), new=TRUE))
  colorbar_horizontal(col=colors_AC, brks=brks, height=4, cex.axis=0.8) # colorbar of potACs in the 4th row
graphics.off()

# Cut of margins:
plot <- image_trim(image_read(plotname))
image_write(plot, path=paste0("cropped_",plotname), format="png")



######### Barplots of PC1 #########

### Constructing indices from 1981 to 2020 just to exemplify creating barplots with base R functions:
pc1_index_ERA <- scale(c(pc1[1:10,1]*(-1), pc1[1:30,3]))[,1]
pc1_index_NCEP1 <- scale(c(pc1[1:10,4], pc1[1:30,6])*(-1))[,1]
pc1_index_NCEP2 <- scale(c(pc1[1:10,7], pc1[1:30,9])*(-1))[,1]
pc_indices <- cbind(pc1_index_ERA, pc1_index_NCEP1, pc1_index_NCEP2) 

### Comparison: correlations and root mean squared errors
rcorr(pc_indices)$r

rmse <- function(o,r) ( sqrt(mean((o-r)**2)) )
for (i in 1:ncol(pc_indices)) {
  print(rmse(pc_indices[,2],pc_indices[,3]))
}

### Plotting data:
# Axis labels:
years <- 1981:2020
years[seq(1,40,2)] <- NA

png(plotname2, units="in", width=30, height=24, res=300, pointsize=64)
par(omi=c(1.3,0,0,0)) # configure outer margin
par(mfrow=c(3,1),
    mar=c(1.8,2.7,1.5,0.2)) # configure margin of each plot
  for (mod in 1:ncol(pc_indices)) {
    barp <- barplot(pc_indices[,mod], ylim=c(-3,3), col=as.character(ifelse(pc1_index>=0, "red", "blue")),
                    xaxs="i", yaxs="i", axes=FALSE, ann=FALSE)
    title(ylab="Standardized index", line=1.8)
    axis(1, at=barp, labels=years, las=2)
    axis(2, at=-3:3, labels=-3:3, las=2)
    text(x=45, y=2.5, label=reans2[mod], adj=1)
  }
graphics.off()

