# calculate the tpi for a set of stations using a 90m DEM.  Buffer should be specified in meters

calculate_station_tpi = function(meta, buffer){

library(raster)
#library(spatialEco)
#library(sp)

m <- raster('Data/DEM/WUS_DEM.tif')
#nlat <- dim(m)[1]
#nlon <- dim(m)[2]
#r<- crop(m, extent(sub_bbox[1],sub_bbox[3],sub_bbox[2],sub_bbox[4]))

# find the grid cell in which each station is located and calculate the tpi at a range of scales:
myPoints <- data.frame(meta$lon, meta$lat)
coordinates(myPoints) <- c('meta.lon', 'meta.lat')
proj4string(myPoints) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

meta$elev - extract(m, myPoints, buffer=buffer, fun=mean, cellnumbers=T) # the buffer units appears to be meters

}


