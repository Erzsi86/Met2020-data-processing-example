####################################################
######### DATA PREPROCESSION & COMPUTATION #########
####################################################

# Monthly geopotential height anomalies are calculated
# in a predefined area and time interval.
# NCEP/NCAR R1 and ERA5 files are handled differently.

corr_pca_monthly_anom <- function(input_path=input_path, input_file=input_file, dataset="ERA5",
                                  lon1=-75, lon2=40, lat1=20, lat2=80, experience_nr=100,
                                  start_date="1981-01-01 00:00:00", end_date="2010-01-31 00:00:00",
                                  output_path=output_path, output_file="absmincor.rds",
                                  output_file2="pca.rds"){
  
 if(dataset=="NCEP") {
    # Load data:
    data <- readRDS(file=paste(input_path, input_file, sep="/"))
   
    # Store dimensions and convert date:
    dim(data) # order of dimensions: latitude, longitude, time
    lat   <- as.numeric(dimnames(data)[[1]])
    lon   <- as.numeric(dimnames(data)[[2]])
    date  <- dimnames(data)[[3]]
    date2 <- as.POSIXct(date, format="%Y_%m_%d_%H", tz="UTC")
  
    # Check dimensions:
    if( length(unique(diff(lon))) == 1 & length(unique(diff(lat))) == 1) {
      print( paste0("Data are available on equidistant grid. Spatial resolution is ",
                    abs(unique(diff(lat))), " x ", abs(unique(diff(lon))), ".") ) } else {
                   "The grid is not equidistant!" }
  
    print(paste0("Data are available between longitudes ", lon[1], " & ", lon[length(lon)],
                 " latitudes ",  lat[1], " & ", lat[length(lat)], ".") )
  
    # Replacing the first two dimensions so the order of the dimensions will be longitude, latitude, time.
    # Then reverting the order of latitudes to get increasing sequence:
    data <- aperm(data, c(2,1,3))
  }
  
  if(dataset=="ERA5") {
    # Load data and store dimensions and convert date:
    nc_file <- nc_open(paste(input_path, input_file, sep="/"))
    lon <- ncvar_get(nc_file, "longitude")
    lat <- ncvar_get(nc_file, "latitude")
    date <- ncvar_get(nc_file, "time")
    date2 <- as.POSIXct(date*3600, origin="1900-01-01 00:00:00.0", tz="UTC")
    data <- ncvar_get(nc_file, "z") / 9.80665 # to get geopotential height [m] rather than geopotential [m^2/s^2]
    
    print(paste0("The examined dataset is: ", input_file))
    
    # Check dimensions:
    if( length(unique(diff(lon))) == 1 & length(unique(diff(lat))) == 1) {
      print( paste0("Data are available on equidistant grid. Spatial resolution is ",
                    abs(unique(diff(lat))), " x ", abs(unique(diff(lon))), ".") ) } else {
                      "The grid is not equidistant!" }
    
    print(paste0("Data are available between longitudes ", lon[1], " & ", lon[length(lon)],
                 " latitudes ",  lat[1], " & ", lat[length(lat)], ".") )
    
    print(paste0("Data are available for ", length(date), " time steps."))
    
    if(diff(lat)[1]<0) {
      data <- data[, length(lat):1, ]
      dims <- dim(data) # order of dimensions: longitude, latitude, time
      lat <- rev(lat)
    }
  }
    
  # Select data:
  lon_sel  <- lon[which(lon==lon1) : which(lon==lon2)]
  lat_sel  <- lat[which(lat==lat1) : which(lat==lat2)]
  start <- which(grepl(start_date, date2))
  end   <- which(grepl(end_date, date2))
  
  data_sel <- data[which(lon==lon1):which(lon==lon2), which(lat==lat1):which(lat==lat2), start:end]
  dims_sel <- dim(data_sel)
  
  # Convert the three-dimensional array to table:
  table <- t(matrix(data=data_sel, nrow=dims_sel[1]*dims_sel[2], ncol=dims_sel[3]))
  
  # Rows indicate times while columns indicate grid cells:
  # in the first column: lat1xlon1, lat1xlon2, ..., lat1xlon144
  # in the last column: lat37xlon1, lat37xlon2, ..., lat25xlon144
  
  # Computing monthly averages:
  group <- format(date2[start:end], format="%Y-%m") # creating grouping variable
  monthly_means <- aggregate(x=table, by=list(group), FUN=mean)
  monthly_means_scaled <- scale(monthly_means[,2:ncol(monthly_means)])
  
  
  #######################################################################
  ######### CORRELATION ANALYSIS & PRINCIPAL COMPONENT ANALYSIS #########
  #######################################################################
  
  ######### COMPUTING STRONGEST NEGATIVE CORRELATIONS #########
  corrs <- rcorr(monthly_means_scaled, type="pearson")

  absmincorrs <- apply(corrs$r, 1, FUN=min)
  absmincorrs_table <- matrix(absmincorrs, nrow=dims_sel[1], ncol=dims_sel[2]) 
  
  saveRDS(object=absmincorrs_table, file=paste(output_path, output_file, sep="/"))

  # Significance test:
  random_corr <- array(NA, dim=c(ncol(monthly_means_scaled),experience_nr))
 
  for (exp in 1:experience_nr) {
    random <- rnorm(nrow(monthly_means_scaled))
    for (col in 1:ncol(monthly_means_scaled)) {
      random_corr[col,exp] <- cor(random, monthly_means_scaled[,col], method="pearson")
    }
  }
  
  random_distribution <- quantile(random_corr, probs=c(0,0.001,0.01,0.05,0.1,0.5,0.9,0.95,0.99,0.999,1))
    
  print(paste0("Distribution of random correlations - number of experiments: ", experience_nr))
  return(round(random_distribution, digits=2))
         
  ######### COMPUTING PRINCIPAL COMPONENT ANALYSIS #########
  res.pca <- prcomp(monthly_means_scaled, center=FALSE, scale=FALSE)
  saveRDS(object=res.pca, file=paste(output_path, output_file2, sep="/"))
}
