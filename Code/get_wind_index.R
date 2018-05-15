
get_wind_index <- function(tab,buffer){

  library(plyr)
  library(dplyr)
  library(ncdf4)
  library(raster)
  library(ctmcmove) # for gradient function
  library(geosphere) # for distance bn coordinates
  library(tidyr)
  
  # set the stations and buffers to use:
  meta <- tab %>% dplyr::select(station_id, lat, lon, elev) %>% distinct() # locations of stations
  #buffer <- c(10000) # in meters, will be the radius around each point
  
  # read in files:
  dem <- raster('Data/DEM/WUS_DEM.tif') # 90m
  dirr <- 'Data/era-i-winds/'
  files <- list.files(dirr, pattern='winds')
  year <- substr(files,6,9)
  
  # set up station points layer to work with dem:
  myPoints <- data.frame(meta$lon, meta$lat)
  coordinates(myPoints) <- c('meta.lon', 'meta.lat')
  proj4string(myPoints) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0")
  
  # set up station points layer to work with winds layer:
  myPointsw <- data.frame(meta$lon+360, meta$lat); names(myPointsw) <- c('meta.lon','meta.lat')
  coordinates(myPointsw) <- c('meta.lon', 'meta.lat')
  proj4string(myPointsw) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0")
  
  # set level to get from 0.75 degree resolution to the buffer resolution:
  # since the distance between lat/lon degrees vary with lat/lon, use the distance in the middle of the stations  
  rlat <- c(mean(meta$lat)-(.75/2), mean(meta$lat)+(.75/2))
  rlon <- c(mean(meta$lon)-(.75/2), mean(meta$lon)+(.75/2))
  
  dlon <- distm(c(rlon[1],rlat[1]), c(rlon[2],rlat[1]), fun = distHaversine) # distance between two adjacent longitudes
  dlat <- distm(c(rlon[1],rlat[1]), c(rlon[1],rlat[2]), fun = distHaversine) # distance between two adjacent latitudes
  #dist_meters <- mean(dlon,dlat)
  
  out <- data.frame(matrix(NA,nrow=(nrow(meta)*length(year)*12),ncol=length(buffer)))
  
  
  for (bb in 1:length(buffer)){
    print(paste0('running buffer ',bb))
  #### GET GRADIENT VECTOR FROM DEM ####
    # aggregate to coarser scale:
    demagg <- aggregate(dem, fact=round(buffer[bb]/90), fun=mean)
    
    # find cell locations of my stations:
    overs <- data.frame(extract(demagg, myPoints, cellnumbers=T))
    
    # get the x and y gradients of the aggregated dem for my point locations:
    grads <- rast.grad(demagg)
    #image(grads$rast.grad.x)
    xgrad <- grads$rast.grad.x[overs$cells]
    ygrad <- grads$rast.grad.y[overs$cells]
  
    gradients <- data.frame(meta,xgrad=xgrad,ygrad=ygrad)
  
  
  #### GET WIND VECTORS ####
    uwindtab <- data.frame(matrix(nrow=nrow(meta), ncol=(length(year)*12)))
    vwindtab <- uwindtab
    for (yy in 1:length(year)){
      print(paste0('year ',yy))
      # open the netcdf file for each year:
      uwind <- brick(paste0(dirr,'winds',year[yy],'.nc'), varname='u') # this coordinate system does match the dem, i checked. also, it is at 0.75 degree res
      vwind <- brick(paste0(dirr,'winds',year[yy],'.nc'), varname='v')
      
    
      # interpolate to a finer grid:
      intlat <- round(dlat/buffer)
      intlon <- round(dlon/buffer)
      #int <- round(dist_meters/buffer)
      #newex <- raster(nrow=dim(uwind)[1]*int, ncol=dim(uwind)[2]*int, ext=extent(uwind))
      newex <- raster(nrow=dim(uwind)[1]*intlat, ncol=dim(uwind)[2]*intlon, ext=extent(uwind))
      
      ufine <- resample(uwind, newex, method='bilinear')
      vfine <- resample(vwind, newex, method='bilinear')
      #image(ufine)
      
      # extract winds for my stations:
      uwindpts <- data.frame(extract(ufine, myPointsw, cellnumbers=T))
      vwindpts <- data.frame(extract(ufine, myPointsw, cellnumbers=T))
      uwindtab[,(1+(yy-1)*12):(yy*12)] <- uwindpts[,2:13]
      names(uwindtab)[(1+(yy-1)*12):(yy*12)] <- names(uwindpts[,2:13])
      vwindtab[,(1+(yy-1)*12):(yy*12)] <- vwindpts[,2:13]
      names(vwindtab)[(1+(yy-1)*12):(yy*12)] <- names(vwindpts[,2:13])
      
    }
    # organize the winds table:
    uwindtab <- cbind(meta, uwindtab)
    vwindtab <- cbind(meta, vwindtab)
    uw <- uwindtab %>% gather('date','uwind', (ncol(meta)+1):ncol(uwindtab))
    vw <- vwindtab %>% gather('date','vwind', (ncol(meta)+1):ncol(vwindtab))
    winds <- left_join(uw,vw, by=c('station_id','lat','lon','elev','date')); rm(uw,vw); gc()
    winds <- winds %>% separate('date', into=c('year','month'), sep=5)
    winds$year <- substr(winds$year,2,5)
    winds$month <- substr(winds$month,2,3)
    
    # combine gradient and winds tables:
    tab <- left_join(gradients, winds, by=c('station_id','lat','lon','elev'))
    tab <- tab %>% mutate(windex = (uwind*xgrad) + (vwind*ygrad)) # yearly and monthly values of windex for each station
    # save all the information, including topographic gradients and wind vectors:
    #write.table(tab, paste0('Data/Windex/windex_',buffer[bb],'_topowx.txt'))
    out[,bb] <- tab$windex
  } # end buffers
  tab <- tab %>% dplyr::select(station_id,year,month)
  out <- bind_cols(tab,out)
  names(out) <- c('station_id','year','month',paste0('windex_',buffer))
  
  # aggregate to mean seasonal values:
  seasons <- data.frame('Fall'=c('09','10','11'), 'Winter'=c('12','01','02'),'Spring'=c('03','04','05'),
                        'Summer'=c('06','07','08'), 'Annual'=c('01','02','03','04','05','06','07','08','09','10','11','12'))
  seasons <- gather(seasons,season,month) %>% distinct()
  out <-out %>% left_join(seasons, by = "month") %>% dplyr::select(-month)
  out <- out %>% group_by(station_id,season,year) %>% summarise_all(mean)
    
  return(out) # returns windex value for each station, year, and season
} # end function


