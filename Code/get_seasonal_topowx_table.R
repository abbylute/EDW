# TopoWx seasonal mean dataset:
library(magrittr)
library(tidyr)
library(plyr)
library(dplyr)
library(tictoc)

# 0. Set these variables:
#.........................................................
datatype <- 'raw' # 'raw','homog',or 'infilled'
varname1 <- paste0('tmax_',datatype)
varname2 <- paste0('tmin_',datatype)
startdate <- '1980-01-01' # these should correspond to calendar years to work with the trend analysis
enddate <- '2016-12-31'
exclude <-  'SNOTEL' #character()# default is
bbox <- c(-125,30,-100,53)
min_pdays <- 0.85  # minimum percent of days of data within each year
min_pyrs <- .90 # minimum percent of years of data at each station for the period defined above
outfn <- paste0('Data/TopoWx/GHCN_',datatype,'_',substr(startdate,1,4),'_summary.txt')


codedir <- 'Code/'
# old way, this was missing the annual months and wasn't an ideal way to do it. actually it accounts for annual separately later on
season_names <- c('Winter','Spring','Summer','Fall')
seasons <- data.frame('Fall'=c('09','10','11'), 'Winter'=c('12','01','02'),'Spring'=c('03','04','05'),
                     'Summer'=c('06','07','08'))
seasons <- gather(seasons,season,month)

#seasons <- c(rep(c('Spring','Summer','Fall','Winter'),each=3),rep('Annual',12))
#seasons <- data.frame(cbind(seasons, c('03','04','05','06','07','08','09','10','11','12','01','02','01','02','03','04','05','06','07','08','09','10','11','12')),stringsAsFactors=F)
#names(seasons) <- c('season','month')

# 1. Get TopoWx data
#.......................................................................
source(paste0(codedir,'get_topowx_data.R'))
tmax <- get_topowx_data(varname=varname1, startdate=startdate, enddate=enddate,
                        bbox=bbox, exclude_network=exclude, min_pdays=min_pdays, min_pyrs=min_pyrs)
tmax_meta <- tmax[[1]]
tmax_meta[,3:5] <- as.numeric(unlist(tmax_meta[,3:5]))
tmax <- tmax[[2]]

tmin <- get_topowx_data(varname=varname2, startdate=startdate, enddate=enddate,
                        bbox=bbox, exclude_network=exclude, min_pdays=min_pdays, min_pyrs=min_pyrs)
tmin_meta <- tmin[[1]]
tmin_meta[,3:5] <- as.numeric(unlist(tmin_meta[,3:5]))
tmin <- tmin[[2]]

# 2. Get seasonal and annual means of Tmax and Tmin datasets
#.......................................................................
cal <- seq(as.Date(startdate),as.Date(enddate),'days')

# TMAX:
  # reorganize the metadata and temperature data into a long data frame:
  temps <- cbind(tmax_meta,tmax)
  tm <- temps %>% gather(key=date, value=temperature, 6:ncol(temps))
  
  datedf <- data.frame(date=c(colnames(temps[6:ncol(temps)])), datenm=cal, stringsAsFactors=F)
  tm <- left_join(tm,datedf, by='date')
  tm <- tm %>% dplyr::select(-date)
  names(tm)[ncol(tm)] <- 'date'
  tm <- tm %>% separate(date, into=c('year','month','day'),sep='-',remove=F)
  tm <- arrange(tm,station_id,date)
  tm <- tm %>% left_join(seasons, by = "month")
  
  # take all the December obs and make year=year+1 then run the season_means so that december gets grouped with winter of the correct year
  tm$year[which(tm$month==12)] <- as.character(as.numeric(tm$year[which(tm$month==12)]) +1)
  tm <- tm %>% filter(!(year==2017 & month==12)) # remove the december 2016 observations which would be grouped with winter 2017 since we don't have all the data for winter 2017
  
  # Calculate seasonal mean temperature values:
  #season_means1 <- tm %>%
  #  group_by(station_id,season,year) %>%
  #  mutate(tmax = mean(temperature,na.rm=T)) %>% dplyr::select(-temperature,-date,-month,-day)
  #season_means1 <- season_means %>% group_by(station_id,season,year) %>% slice(which.min(tmax))
  
  season_means <- tm %>%
    group_by(station_id,season,year,network,lat,lon,elev) %>%
    summarise(tmax = case_when(
      length(which(is.na(temperature))) > 20 ~ as.numeric(NA), # if a given season is missing more than 20 daily values then don't calculate a seasonal mean it
      length(which(is.na(temperature))) <= 20 ~ mean(temperature,na.rm=T))) 
  

  # Calculate water year mean temperature values:
  # first, bump oct and nov to the next year (already did this for december above)
  tma <- tm
  tma$year[which(tma$month %in% c(10,11))] <- as.character(as.numeric(tma$year[which(tma$month %in% c(10,11))]) +1)
  tma <- tma %>% filter(!(year==2017 & month %in% c(10,11))) # remove the oct and nov 2016 observations which would be grouped with wy 2017 since we don't have all the data for wy 2017
  
  #annual_means1 <- tma %>% group_by(station_id,year) %>% mutate(tmax=mean(temperature, na.rm=T)) %>%
  #  mutate(season='Annual') %>% dplyr::select(-c(month,day,date,temperature)) %>% distinct()
 
  annual_means <- tma %>% group_by(station_id,year,network,lat,lon,elev) %>% 
    summarise(tmax = case_when(
      length(which(is.na(temperature))) > 60 ~ as.numeric(NA),
      length(which(is.na(temperature))) <= 60 ~ mean(temperature, na.rm=T))) %>%
    mutate(season='Annual') 
  
  tmax_master <- full_join(season_means, annual_means, by=c('station_id','network','lat','lon','elev','year','season','tmax'))
  rm(tm,tma,season_means,annual_means,tmax,tmax_meta); gc()
  
  # TMIN
  # reorganize the metadata and temperature data into a long data frame:
  temps <- cbind(tmin_meta,tmin)
  tm <- temps %>% gather(key=date, value=temperature, 6:ncol(temps))
  
  datedf <- data.frame(date=c(colnames(temps[6:ncol(temps)])), datenm=cal, stringsAsFactors=F); rm(temps); gc()
  tm <- left_join(tm,datedf, by='date')
  tm <- tm %>% dplyr::select(-date)
  names(tm)[ncol(tm)] <- 'date'
  tm <- tm %>% separate(date, into=c('year','month','day'),sep='-',remove=F)
  tm <- arrange(tm,station_id,date)
  tm <- tm %>% left_join(seasons, by = "month")
  
  # take all the December obs and make year=year+1 then run the season_means so that december gets grouped with winter of the correct year
  tm$year[which(tm$month==12)] <- as.character(as.numeric(tm$year[which(tm$month==12)]) +1)
  tm <- tm %>% filter(!(year==2017 & month==12)) # remove the december 2016 observations which would be grouped with winter 2017 since we don't have all the data for winter 2017
  
  # Calculate seasonal mean temperature values:
  #season_means <- tm %>%
  #  group_by(station_id,season,year) %>%
  #  mutate(tmin = mean(temperature,na.rm=T)) %>% dplyr::select(-temperature,-date,-month,-day)
  #season_means <- season_means %>% group_by(station_id,season,year) %>% slice(which.min(tmin))
  
  season_means <- tm %>%
    group_by(station_id,season,year,network,lat,lon,elev) %>%
    summarise(tmin = case_when(
      length(which(is.na(temperature))) > 20 ~ as.numeric(NA), # if a given season is missing more than 20 daily values then don't calculate a seasonal mean it
      length(which(is.na(temperature))) <= 20 ~ mean(temperature,na.rm=T))) 
  
  # Calculate water year mean temperature values:
  # first, bump oct and nov to the next year (already did this for december above)
  tma <- tm
  tma$year[which(tma$month %in% c(10,11))] <- as.character(as.numeric(tma$year[which(tma$month %in% c(10,11))]) +1)
  tma <- tma %>% filter(!(year==2017 & month %in% c(10,11))) # remove the oct and nov 2016 observations which would be grouped with wy 2017 since we don't have all the data for wy 2017
  
  #annual_means <- tma %>% group_by(station_id,year) %>% mutate(tmin=mean(temperature, na.rm=T)) %>%
  #  mutate(season='Annual') %>% dplyr::select(-c(month,day,date,temperature)) %>% distinct()
  
  annual_means <- tma %>% group_by(station_id,year,network,lat,lon,elev) %>% 
    summarise(tmin = case_when(
      length(which(is.na(temperature))) > 60 ~ as.numeric(NA),
      length(which(is.na(temperature))) <= 60 ~ mean(temperature, na.rm=T))) %>%
    mutate(season='Annual') 
  
  tmin_master <- full_join(season_means, annual_means, by=c('station_id','network','lat','lon','elev','year','season','tmin'))
  rm(tm,tma,season_means,annual_means,tmin,tmin_meta); gc()
  
  
# 3. Join the Tmax and Tmin tables:
  master <- full_join(tmax_master,tmin_master, by=c('station_id','network','lat','lon','elev','year','season'))
  rm(tmax_master,tmin_master,get_topowx_data);gc()
  
  ######## new code:
  #source('Code/prep_station_data.R')
  #dat <- prep_station_data(master, save_dir='Data/TopoWx/GHCN_raw_1980_summary_wi.txt')
    
  # actually, just run get_wind_index for now
  #source('Code/get_wind_index.R')
  buffer <- c(500,1000,5000,10000,20000,40000)
  for (bb in 1:length(buffer)){
    tic(paste0('buffer ',buffer[bb]))
    out <- get_wind_index(master,buffer[bb],save_dir='Data/Windex/')
    write.table(out, paste0('Data/Windex/windex_',buffer[bb],'_topowx_GHCN_raw.txt'))
    toc()
  }
  
  
  ######## old code:
# 4. Add TPI data:
  trim_data <- master %>% group_by(station_id,add=F)%>% select(station_id,lat,lon,elev)  %>% distinct()
  
  source('Code/calculate_station_tpi.R')
  buffers <- c(200,500,1000,5000)
  tic('calculate tpis') # 30 minutes
  tpis <- calculate_station_tpi(trim_data,buffers)
  toc()
  tpis <- cbind(trim_data$station_id,tpis); names(tpis)[1] <- 'station_id'
  master <- left_join(master,tpis,by='station_id')
  
  
# 5. Add Solar Radiation Data:
  source('Code/get_station_srad.R')
  tic('get solar')
  srad <- get_station_srad(trim_data) 
  toc()
  srad <- cbind(trim_data$station_id,srad); names(srad)[1] <- 'station_id'
  
  # convert monthly srad to seasonal srad:
  srad_seasonal <- srad %>% group_by(station_id) %>% 
    summarise(Winter = mean(c(jan,feb,dec),na.rm=T), Spring = mean(c(mar,apr,may),na.rm=T),
              Summer = mean(c(jun,jul,aug),na.rm=T), Fall = mean(c(sep,oct,nov),na.rm=T),
              Annual = mean(c(jan,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec))) %>%
    gather('season','srad',c('Winter','Spring','Summer','Fall','Annual'))
  
  master <- left_join(master,srad_seasonal, by=c('station_id','season'))
  
# 6. Save table
  write.table(master, outfn)

