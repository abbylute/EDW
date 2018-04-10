# get solar radiation values for station locations
# this grabs topographically corrected srad data that John created (04/2018) using my 90m DEM

# returns a dataframe with each row corresponding to a station in meta and each column corresponding to a monthly srad value

get_station_srad = function(meta){
  
  if (!all(c('lat','lon','elev') %in% names(meta))){
    print('dataframe called meta must contain column names lat, lon, and elev')
  }
  
  library(ncdf4)
  dirr <- 'Data/Solar_radiation/'
  
  # preallocate the output dataframe:
  srad_out <- data.frame(matrix(nrow=nrow(meta), ncol=12))
  names(srad_out) <- c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
  
  # extract srad for each month:
  for (mm in 1:12){ 
    nc <- nc_open(paste0(dirr,'srad_',mm,'.nc'), write=F)
    lats <- ncvar_get(nc,'lat')
    lons <- ncvar_get(nc,'lon')
    
    plats <- numeric()
    plons <- numeric()
    for (ss in 1:nrow(meta)){
      plats[ss] <- which.min(abs(lats-meta$lat[ss]))
      plons[ss] <- which.min(abs(lons-meta$lon[ss]))
    }
    
    # just extract srad values starting with the lowest lat and lon we need so we don't have to grab the entire thing
    srad <- ncvar_get(nc,'srad', start=c(min(plons),min(plats)), count=c(-1,-1)) # [lon,lat]
    
    # reset the indices to extract since above you extracted starting from first station index:
    plons <- plons-min(plons)+1
    plats <- plats-min(plats)+1
    
    # extract the station srad values:
    srad_out[,mm] <- srad[cbind(plons,plats)]
    
  }
  
  return(srad_out)
}



