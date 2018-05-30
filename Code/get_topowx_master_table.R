# This script calculates yearly and seasonal values of tmax and tmin for every station in the TopoWx dataset.  
# It also calculates all of the TPI, solar radiation, and windices for these stations.  At the moment I've 
# only run this over the raw temperature data.  To do the same for the homogenized or infilled temperature 
# data, I would only need to re-run the temperature part, not the indices.

library(ncdf4)
library(tidyverse)
library(tictoc)
library(parallel)

# this accounts for annual separately later on
season_names <- c('Winter','Spring','Summer','Fall')
seasons <- data.frame('Fall'=c('09','10','11'), 'Winter'=c('12','01','02'),'Spring'=c('03','04','05'),
                      'Summer'=c('06','07','08'))
seasons <- gather(seasons,season,month)


# Get all topowx tmax stations:
  ncfile <- 'Data/TopoWx/stn_obs_tmax.nc'
  nc <- nc_open(ncfile,write=F)
  lat <- ncvar_get(nc, 'latitude')
  lon <- ncvar_get(nc, 'longitude')
  station_id <- ncvar_get(nc, 'station_id')
  elev <- ncvar_get(nc, 'elevation')
  tmax_meta <- bind_cols(station_id=station_id,lat=lat,lon=lon,elev=elev)
  
  cal <- seq(as.Date('1948-01-01'),as.Date('2016-12-31'),'days')
  tmax_raw <- ncvar_get(nc, 'tmax_raw', start=c(1,1), count=c(-1,-1))
  tmax_raw <- cbind(tmax_meta,tmax_raw); names(tmax_raw) <- c('station_id','lat','lon','elev',as.character(cal))
  #tmax_hom <- ncvar_get(nc, 'tmax_homog', start=c(1,1), count=c(-1,-1))
  #tmax_inf <- ncvar_get(nc, 'tmax_infilled', start=c(1,1), count=c(-1,-1))
  
  tm <- tmax_raw %>% gather(key=date, value=temperature, 5:ncol(tmax_raw))
  
  tm$year <- substr(tm$date,1,4)
  tm$month <- substr(tm$date,6,7)
  tm$day <- substr(tm$date, 9,10)
  tm <- arrange(tm,station_id,date)
  tm <- tm %>% left_join(seasons, by = "month")
  
  # take all the December obs and make year=year+1 then run the season_means so that december gets grouped with winter of the correct year
  tm$year[which(tm$month==12)] <- as.character(as.numeric(tm$year[which(tm$month==12)]) +1)
  tm <- tm %>% filter(!(year==2017 & month==12)) # remove the december 2016 observations which would be grouped with winter 2017 since we don't have all the data for winter 2017
  
  season_means <- tm %>%
    group_by(station_id,season,year,lat,lon,elev) %>%
    summarise(tmax = case_when(
      length(which(is.na(temperature))) > 20 ~ as.numeric(NA), # if a given season is missing more than 20 daily values then don't calculate a seasonal mean it
      length(which(is.na(temperature))) <= 20 ~ mean(temperature,na.rm=T))) 
  
  # Calculate water year mean temperature values:
  # first, bump oct and nov to the next year (already did this for december above)
  tm$year[which(tm$month %in% c(10,11))] <- as.character(as.numeric(tm$year[which(tm$month %in% c(10,11))]) +1)
  tm <- tm %>% filter(!(year==2017 & month %in% c(10,11))) # remove the oct and nov 2016 observations which would be grouped with wy 2017 since we don't have all the data for wy 2017
  
  annual_means <- tm %>% group_by(station_id,year,lat,lon,elev) %>% 
    summarise(tmax = case_when(
      length(which(is.na(temperature))) > 60 ~ as.numeric(NA),
      length(which(is.na(temperature))) <= 60 ~ mean(temperature, na.rm=T))) %>%
    mutate(season='Annual') 
  
  tmax_master <- full_join(season_means, annual_means, by=c('station_id','lat','lon','elev','year','season','tmax'))
  write.table(tmax_master,'Data/TopoWx/All_topowx_tmax.txt')
  rm(tm,season_means,annual_means,tmax_raw,tmax_meta); gc()
 
  
# Get all topowx tmin stations:
  ncfile <- 'Data/TopoWx/stn_obs_tmin.nc'
  nc <- nc_open(ncfile,write=F)
  lat <- ncvar_get(nc, 'latitude')
  lon <- ncvar_get(nc, 'longitude')
  station_id <- ncvar_get(nc, 'station_id')
  elev <- ncvar_get(nc, 'elevation')
  tmin_meta <- bind_cols(station_id=station_id,lat=lat,lon=lon,elev=elev)
  
  cal <- seq(as.Date('1948-01-01'),as.Date('2016-12-31'),'days')
  tmin_raw <- ncvar_get(nc, 'tmin_raw', start=c(1,1), count=c(-1,-1))
  tmin_raw <- cbind(tmin_meta,tmin_raw); names(tmin_raw) <- c('station_id','lat','lon','elev',as.character(cal))
  #tmax_hom <- ncvar_get(nc, 'tmax_homog', start=c(1,1), count=c(-1,-1))
  #tmax_inf <- ncvar_get(nc, 'tmax_infilled', start=c(1,1), count=c(-1,-1))
  
  tm <- tmin_raw %>% gather(key=date, value=temperature, 5:ncol(tmin_raw))
  
  tm$year <- substr(tm$date,1,4)
  tm$month <- substr(tm$date,6,7)
  tm$day <- substr(tm$date, 9,10)
  tm <- arrange(tm,station_id,date)
  tm <- tm %>% left_join(seasons, by = "month")
  
  # take all the December obs and make year=year+1 then run the season_means so that december gets grouped with winter of the correct year
  tm$year[which(tm$month==12)] <- as.character(as.numeric(tm$year[which(tm$month==12)]) +1)
  tm <- tm %>% filter(!(year==2017 & month==12)) # remove the december 2016 observations which would be grouped with winter 2017 since we don't have all the data for winter 2017
  
  season_means <- tm %>%
    group_by(station_id,season,year,lat,lon,elev) %>%
    summarise(tmin = case_when(
      length(which(is.na(temperature))) > 20 ~ as.numeric(NA), # if a given season is missing more than 20 daily values then don't calculate a seasonal mean it
      length(which(is.na(temperature))) <= 20 ~ mean(temperature,na.rm=T))) 
  
  # Calculate water year mean temperature values:
  # first, bump oct and nov to the next year (already did this for december above)
  tm$year[which(tm$month %in% c(10,11))] <- as.character(as.numeric(tm$year[which(tm$month %in% c(10,11))]) +1)
  tm <- tm %>% filter(!(year==2017 & month %in% c(10,11))) # remove the oct and nov 2016 observations which would be grouped with wy 2017 since we don't have all the data for wy 2017
  
  annual_means <- tm %>% group_by(station_id,year,lat,lon,elev) %>% 
    summarise(tmin = case_when(
      length(which(is.na(temperature))) > 60 ~ as.numeric(NA),
      length(which(is.na(temperature))) <= 60 ~ mean(temperature, na.rm=T))) %>%
    mutate(season='Annual') 
  
  tmin_master <- full_join(season_means, annual_means, by=c('station_id','lat','lon','elev','year','season','tmin'))
  write.table(tmin_master,'Data/TopoWx/All_topowx_tmin.txt')
  rm(tm,season_means,annual_means,tmin,tmin_meta); gc()
  
  nc_close(nc); gc()
  
# join tmax and tmin tables
  topowxtab <- full_join(tmax_master,tmin_master, by = c('station_id','lat','lon','elev','year','season'))
  class(topowxtab$lat) <- 'numeric'
  class(topowxtab$lon) <- 'numeric'
  class(topowxtab$elev) <- 'numeric'
  
  topowxtab <- topowxtab %>% filter(lon < -100) # only WUS

  rm(tmaxtab,tmintab,elev,lat,lon,station_id,nc,ncfile); gc()
  
# Add TPI data:

  source('Code/calculate_station_tpi.R')
  buffers <- c(200,500,1000,5000)
  tic('calculate tpis') # took an 1.5 hours
  tpis <- calculate_station_tpi(topowxtab,buffers)
  toc()
  
  master <- cbind(topowxtab,tpis)
  
  
# Add Solar Radiation Data:
  
  source('Code/get_station_srad.R')
  tic('get solar') # about 30 minutes
  srad <- get_station_srad(topowxtab) # outputs seasonal values of srad_cl and srad_po
  toc()
  
  master <- left_join(master,srad, by=c('station_id','lat','lon','season'))
  write.table(master,'Data/TopoWx/All_topowx_table_part1.txt')
  
# Get Windices:  
  # this first loop took almost 5 hours
  source('Code/get_wind_index.R')
  buffer <- rev(c(500,1000,5000,10000,20000,40000))
  for (bb in 1:length(buffer)){
    tic(paste0('buffer ',buffer[bb]))
    out <- get_wind_index(master,buffer[bb],save_dir='Data/Windex/')
    write.table(out, paste0('Data/Windex/windex_',buffer[bb],'_all_topowx.txt'))
    toc()
  }
  
  # join all the wind files together:
  for (bb in 1:length(buffer)){
    if (bb==1){
      wind <- read.table(paste0('Data/Windex/windex_',buffer[bb],'_all_topowx.txt'))
    } else {
      windadd <- read.table(paste0('Data/Windex/windex_',buffer[bb],'_all_topowx.txt'))
      wind <- left_join(wind,windadd, by= c('station_id','season','year'))
    }
  }
  wind$year <- as.character(wind$year)
  
  master <- left_join(master,wind, by=c('station_id','season','year'))
  #write.table(master, 'Data/TopoWx/All_topowx_table_part2')
  
  
# Get free-air correlations
  # skip this for now since I'm not sure exactly how to do it
  #source('Code/get_free_air_temperatures.R')
  
  #trim_master <- master %>% dplyr::select(station_id,lat,lon) %>% distinct()
  #level <- '700'
  #radius <- 1 # in degrees
  
  # for each station:
  #st_bbox <- cbind(trim_master$lon-radius, trim_master$lat-radius, trim_master$lon+radius, trim_master$lat+radius)
  
  #cl <- makeCluster(7)
  #clusterExport(cl, c("get_free_air_temperatures","st_bbox","level"))
  #tic('go')
  #fa <- parApply(cl,st_bbox, 1, function(i) get_free_air_temperatures(bbox=i,level))
  #toc()
  #stopCluster(cl)
  
  # reorganize the list of dfs into a large df
  #fadf <- dplyr::bind_rows(fa)
  
  #trim_master_long <- bind_cols(station_id=rep(trim_master$station_id, each=nrow(fa[[1]])))
  #fadf <- bind_cols(trim_master_long, fadf)
  #fadf$year <- as.integer(fadf$year)
  
  # Compute correlations between free-air and station temperatures across years:
  #fa_master <- left_join(master, fadf, by =c('station_id','year','season'))
  #fa_master <- fa_master %>% filter(!is.na(fatemp))
  #fa_out <- fa_master %>% group_by(station_id,season) %>% 
  #  summarise(fa_tmax_cor = cor.test(tmax,fatemp)$estimate,
  #            fa_tmin_cor = cor.test(tmin,fatemp)$estimate)
  
  #master <- left_join(master, fa_out, by=c('station_id','season'))
  
  
  
# Save table:
  write.table(master, 'Data/TopoWx/All_topowx_master.txt')
  save(master, file='Data/TopoWx/All_topowx_master.RData')
  
  