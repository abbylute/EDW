# brks should include outside extents


get_elev_transects = function(brks, orient = 'horizontal', fun=mean, bbox = c(-125,30,-100,53)){
  library(raster)
  brks <- sort(brks)
  midbrks <- (brks[1:(length(brks)-1)] + brks[2:length(brks)])/2
  
  m <- raster('Documents/IGERT/Data/DEM/WUS_DEM.tif')
  nlat <- dim(m)[1]
  nlon <- dim(m)[2]

  out <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("lat", "lon", "elev"))
  
  if (orient== 'horizontal'){
    for (bb in 1:(length(brks)-1)){
      # extract the band:
      r <- crop(m, extent(bbox[1],bbox[3],brks[bb],brks[bb+1]))
      # aggregate over columns (lats)
      aa <- data.frame(rasterToPoints(aggregate(r, fact=c(1,dim(r)[1]) , fun=fun, na.rm=T)))
      out <- rbind(out, aa)
    }
    names(out)<- c('lon','lat','elev')
    out$facet = factor(out$lat, levels = sort(unique(out$lat), decreasing=T), labels=rev(midbrks))
  } else {  # if looking at vertical  (longitude bands)
    for (bb in 1:(length(brks)-1)){
      # extract the band:
      r <- crop(m, extent(brks[1],brks[2],bbox[2],bbox[4]))
      # aggregate over columns (lats)
      aa <- data.frame(rasterToPoints(aggregate(r, fact=c(dim(r)[2], 1) , fun=fun, na.rm=T)))
      out <- rbind(out, aa)
    } 
    names(out)<- c('lon','lat','elev')
    out$facet = factor(out$lat, levels = sort(unique(out$lat), decreasing=F), labels=midbrks)
  }
  return(out)
}

