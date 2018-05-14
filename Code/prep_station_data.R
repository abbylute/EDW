# Prepare Station Data:

# takes a dataframe of station metadata and grabs the tpi and srad values for those stations

# data should contain columns for station_id, lat, lon, elev, and month (either 'jan' or '01')

# save_dir can be either F (default) or a file name to save the resulting data to.

prep_station_data = function(data, save_dir=F){
  
  if (!all(c('station_id','lat','lon','elev','month') %in% names(data))){
    print('Error: dataframe called data must contain columns station_id, lat, lon, elev, and month')
  }
  
  trim_data <- data %>% dplyr::select(station_id,lat,lon,elev) %>% distinct()
  
  print('Grabbing TPI values...')
  
  # Get TPI values:
  source('Code/calculate_station_tpi.R')
  buffers <- c(200,500,1000,5000)
  tpis <- calculate_station_tpi(trim_data,buffers)
  tpis <- cbind(trim_data$station_id,tpis); names(tpis)[1] <- 'station_id'
  data <- left_join(data,tpis,by='station_id')
  
  print('Grabbing seasons...')
  # Get seasons:
  if (is.numeric(data$month) | nchar(data$month[1])==2){
    seasons <- data.frame('Fall'=c('09','10','11'), 'Winter'=c('12','01','02'),'Spring'=c('03','04','05'),
                          'Summer'=c('06','07','08'), , 'Annual'=c('01','02','03','04','05','06','07','08','09','10','11','12'))
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
  # take all the December obs and make year=year+1 then run the season_means so that december gets grouped with winter of the correct year
  #tm$year[which(tm$month==12)] <- as.character(as.numeric(tm$year[which(tm$month==12)]) +1)
  #tm <- tm %>% filter(!(year==2017 & month==12)) # remove the december 2016 observations which would be grouped with winter 2017 since we don't have all the data for winter 2017
  
  print('Grabbing solar radiation values...')
  # Get solar radiation:
  source('Code/get_station_srad.R')
  srad <- get_station_srad(trim_data) 
  srad <- cbind(trim_data$station_id,srad); names(srad)[1] <- 'station_id'
  
  # convert monthly srad to seasonal srad:
  srad_seasonal <- srad %>% group_by(station_id) %>% 
    summarise(Winter = mean(c(jan,feb,dec),na.rm=T), Spring = mean(c(mar,apr,may),na.rm=T),
              Summer = mean(c(jun,jul,aug),na.rm=T), Fall = mean(c(sep,oct,nov),na.rm=T)) %>%
    gather('season','srad',c('Winter','Spring','Summer','Fall'))
  
  if ('srad' %in% names(data)){ # in case there is an existing 'srad' column in the dataframe
    names(data)[which(names(data)=='srad')] <- 'srad_orig'
  }
  
  data <- left_join(data,srad_seasonal, by=c('station_id','season'))
  
  print('Grabbing windices...')
  # Get windices:
  source('Code/get_wind_index.R')
  buffers <- c(500, 1000, 5000, 10000, 30000) # in meters
  wi <- get_wind_index(data,buffers)
  class(wi$year) <- 'integer'
  data <- left_join(data,wi,by=c('station_id','season','year'))
  
  
  
  
  if (!save_dir == F){
    write.table(data, save_dir)
  }
  
  return(data)
  
}