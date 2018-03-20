# calculate the tpi for a set of stations using a 90m DEM.  Buffer should be specified in meters

calculate_station_tpi = function(meta, buffer){

  library(raster)
  
  m <- raster('Data/DEM/WUS_DEM.tif')
  
  # find the grid cell in which each station is located and calculate the tpi at a range of scales:
  myPoints <- data.frame(meta$lon, meta$lat)
  coordinates(myPoints) <- c('meta.lon', 'meta.lat')
  proj4string(myPoints) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  
  out <- data.frame(matrix(NA,nrow=nrow(meta),ncol=length(buffer)))

  for (ii in 1:length(buffer)){
    out[,ii] <- meta$elev - extract(m, myPoints, buffer=buffer[ii], fun=mean, cellnumbers=T) # the buffer units appears to be meters
  }
    names(out) <- paste0('m',buffer)
    out
}


