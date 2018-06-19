# Prepare Station Data:

# takes a dataframe of station metadata and grabs the tpi and srad values for those stations

# data should contain columns for station_id, lat, lon, elev, and month (either 'jan' or '01')

# save_dir can be either F (default) or a file name to save the resulting data to.

prep_station_data = function(data, save_dir=F){
  
  if (!all(c('station_id','lat','lon','elev') %in% names(data))){
    print('Error: dataframe called data must contain columns station_id, lat, lon, elev, and month')
  }
  if (! ('season' %in% names(data) | 'month' %in% names(data))){
    print('Error: dataframe called data must either contain a month column or a season column')
  }
  
  library(raster)
  
  
  trim_data <- data %>% group_by(add=F) %>% dplyr::select(station_id,lat,lon,elev)  %>% distinct()
  class(trim_data$lat) <- 'numeric'
  class(trim_data$lon) <- 'numeric'
  class(trim_data$elev) <- 'numeric'
  
  
  # Get TPI values:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  print('Grabbing TPI values...')
  
  source('Code/calculate_station_tpi.R')
  buffers <- c(200,500,1000,5000)
  tpis <- calculate_station_tpi(trim_data,buffers)
  tpis <- cbind(trim_data$station_id,tpis,stringsAsFactors=F); names(tpis)[1] <- 'station_id'
  data <- left_join(data,tpis,by='station_id')
  rm(tpis,buffers); gc()
  
  
  # Get Aspect and Slope::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  print('Grabbing Aspect and Slope values...')
  
  rast <- brick('Data/DEM/WUS_DEM_slope_aspect.tif')
  
  myPoints <- data.frame(trim_data$lon, trim_data$lat)
  coordinates(myPoints) <- c('trim_data.lon', 'trim_data.lat')
  proj4string(myPoints) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") # proj4string(rast)
  
  out <- as.data.frame(raster::extract(rast, myPoints)) 
  out <- bind_cols('station_id' = trim_data$station_id,out); names(out) <- c('station_id','slope','aspect')
  data <- left_join(data,out,by='station_id')
  
  rm(rast,out); gc()
  
  
  # Get seasons:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  print('Grabbing seasons...')

  if ('season' %in% names(data)){
    print('data already has seasons')
  } else {
    if (is.numeric(data$month) | nchar(data$month[1])==2){
      seasons <- data.frame('Fall'=c('09','10','11'), 'Winter'=c('12','01','02'),'Spring'=c('03','04','05'),
                            'Summer'=c('06','07','08'), 'Annual'=c('01','02','03','04','05','06','07','08','09','10','11','12'))
      seasons <- gather(seasons,season,month) %>% distinct()
      data <- data %>% left_join(seasons, by = "month")
    } else if (nchar(data$month[1])==3) {
      seasons <- data.frame('Fall'=c('sep','oct','nov'), 'Winter'=c('dec','jan','feb'),'Spring'=c('mar','apr','may'),
                            'Summer'=c('jun','jul','aug'), 'Annual'=c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'))
      seasons <- gather(seasons,season,month) %>% distinct()
      data <- data %>% left_join(seasons, by = "month")
    } else {
      print('Error: not sure how to convert month field to seasons.  Try one of these formats: jan or 01')
    }
  }
  # take all the December obs and make year=year+1 then run the season_means so that december gets grouped with winter of the correct year
  #tm$year[which(tm$month==12)] <- as.character(as.numeric(tm$year[which(tm$month==12)]) +1)
  #tm <- tm %>% filter(!(year==2017 & month==12)) # remove the december 2016 observations which would be grouped with winter 2017 since we don't have all the data for winter 2017
  
  
  # Get solar radiation::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  print('Grabbing solar radiation values...')

  source('Code/get_station_srad.R')
  srad <- get_station_srad(trim_data) 

  if ('srad' %in% names(data)){ # in case there is an existing 'srad' column in the dataframe
    names(data)[which(names(data)=='srad')] <- 'srad_orig'
  }
  
  data <- left_join(data,srad, by=c('station_id','lat','lon','season'))
  rm(srad); gc()
  
  
  # Get windices::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  print('Grabbing windices...')

  #source('Code/get_wind_index.R')
  winddir <- 'Data/Windex/'
  buffer <- c(500,1000,5000,10000,20000,40000)# in meters
  tic('ready go')
  for (bb in 1:length(buffer)){
    print(paste0('running ',buffer[bb],'m buffer'))
    tic('start buffer')
    wind_sub <- read.table(paste0(winddir,'windex_',buffer[bb],'.txt'), stringsAsFactors=F)
    wind_sub <- wind_sub %>% filter(station_id %in% trim_data$station_id)
    if (!bb == 1){
      winds <- left_join(winds,wind_sub, by=c('station_id','season','year'))
    } else {
      winds <- wind_sub
    }
    toc()
  }
  toc()
  
  data <- left_join(data,winds,by=c('station_id','season','year'))
  rm(winddir,buffer,wind_sub,winds,bb); gc()
  
  
  # Get cloud cover::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  print('Grabbing cloud cover...')
  cdir <- 'Data/MODIS_clouds/'
  
  seas <- unique(data$season)
  clouds <- data.frame(matrix(nrow=nrow(trim_data), ncol=length(seas))); names(clouds) <- seas
  for (ss in 1:length(seas)){
    print(seas[ss])
    rast <- raster(paste0(cdir,tolower(seas[ss]),'_cloudcover.tif'))
    clouds[,ss] <- raster::extract(rast, myPoints)
  }
  clouds <- bind_cols(trim_data,clouds)
  clouds <- clouds %>% gather('season','cloudcover',seas)
  
  data <- left_join(data,clouds, by=c('station_id','lat','lon','elev','season'))
  rm(cdir,clouds,rast,ss); gc()
  
  
  # Get land cover::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  print('Grabbing land cover...')
  
  rast <- brick('Data/EarthEnv_landcover/WUS_landcover.tif')
  landcover <- as.data.frame(raster::extract(rast, myPoints)) 
  names(landcover) <- paste0('landcover_',c(1:12))
  landcover <- bind_cols(trim_data, landcover)              
  
  data <- left_join(data, landcover, by=c('station_id','lat','lon','elev'))
  rm(rast,landcover); gc()
  
  
  # Get waterbody distance index ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  print('Grabbing water body distance index...')
  source('Code/prep_waterbody_layer.R')
  
  water <- prep_waterbody_layer(trim_data)
  
  data <- left_join(data,water, by=c('station_id','lat','lon','elev'))
  rm(water); gc()
  
  
  # Get free-air temperature ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  print('Grabbing free-air temperature...')
  
  fatemp <- get(load('Data/era-i-heights/free_air_temperature_yearly_seasonally.RData'))
  # note the above file was created in the script 'get_free_air_adjusted_station_temp...'
  fa_lons <- unique(fatemp$lon)
  fa_lats <- unique(fatemp$lat)
  
  stations <- trim_data %>% group_by(station_id,lat,lon) %>% mutate(fa_lon = fa_lons[which.min(abs(fa_lons-lon))],
                                                                   fa_lat = fa_lats[which.min(abs(fa_lats-lat))])
  stations <- stations %>% group_by(add=F) %>% select(station_id,fa_lon,fa_lat)
  fa_master <- left_join(stations,fatemp,by=c("fa_lat"="lat", "fa_lon"="lon")) %>% select(-fa_lon,-fa_lat)
  class(fa_master$station_id) <- 'character'
  class(fa_master$year) <- 'integer'
  class(fa_master$season) <- 'character'
  
  data <- left_join(data,fa_master, by=c('station_id','year','season'))
  data$fa_adj_temp <- data$temp-data$fatemp
  rm(fatemp,fa_lons,fa_lats,stations,fa_master); gc()
    
    
    
  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  if (!save_dir == F){
    write.table(data, save_dir)
  }
  
  return(data)
  
}