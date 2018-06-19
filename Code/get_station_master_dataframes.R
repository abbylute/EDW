# Prepare station dataset for analysis:
# trim to complete records
# get all indices


startdate <- '1979-10-01'
enddate <- '2016-09-30'
min_pdays <- 0.85
min_pyrs <- 0.90
bbox <- c(-125, 30, -100, 53)
outdir <- 'Data/Station_dataframes/'

cal <- seq(as.Date(startdate),as.Date(enddate),'days')
source('Code/get_topowx_data.R')
source('Code/topowx_data_to_seasonal.R')
source('Code/ghcn_data_to_complete_seasonal.R')
source('Code/prep_station_data.R')
library(tidyr)
library(plyr)
library(dplyr)
library(tictoc)

# Grab datasets:
#------------------------------------------------------------------
# get topowx datasets:
  # tmax
  tmax_topo <- get_topowx_data('tmax_homog',startdate,enddate,bbox,exclude_network=character(),min_pdays,min_pyrs)
  tmax_topo_meta <- tmax_topo[[1]] %>% select(-network)
  tmax_topo <- tmax_topo[[2]]
  
  # convert to seasonal values:
  tmax_topo_seas <- topowx_data_to_seasonal(tmax_topo_meta,tmax_topo,cal)
  
  #save this intermediate step:
  write.table(tmax_topo_seas,paste0(outdir,'Topowx_homog_complete_seasonal_tmax.txt'))

  # tmin
  tmin_topo <- get_topowx_data('tmin_homog',startdate,enddate,bbox,exclude_network=character(),min_pdays,min_pyrs)
  tmin_topo_meta <- tmin_topo[[1]] %>% select(-network)
  tmin_topo <- tmin_topo[[2]]
  
  # convert to seasonal values:
  tmin_topo_seas <- topowx_data_to_seasonal(tmin_topo_meta,tmin_topo,cal)
  
  #save this intermediate step:
  write.table(tmin_topo_seas,paste0(outdir,'Topowx_homog_complete_seasonal_tmin.txt'))

  rm(tmax_topo,tmax_topo_meta,tmin_topo,tmin_topo_meta); gc()


# get GHCND datasets:
  # tmax
  tmax_ghcn <- read.table('Data/GHCN_daily/Tmax_WUS_ghcnd_clean.txt')
  
  # get seasonal and complete data:
  tic('get ghcn seasonal') # about 16 minutes
  tmax_ghcn_seas <- ghcn_data_to_complete_seasonal(tmax_ghcn,cal,min_pdays,min_pyrs)
  toc()
  
  #save this intermediate step:
  write.table(tmax_ghcn_seas,paste0(outdir,'GHCND_clean_complete_seasonal_tmax.txt'))
  rm(tmax_ghcn); gc()
  
  # tmin
  tmin_ghcn <- read.table('Data/GHCN_daily/Tmin_WUS_ghcnd_clean.txt')
  
  # get seasonal and complete data:
  tic('get ghcn seasonal') # about 14 minutes
  tmin_ghcn_seas <- ghcn_data_to_complete_seasonal(tmin_ghcn,cal,min_pdays,min_pyrs)
  toc()
  
  #save this intermediate step:
  write.table(tmin_ghcn_seas,paste0(outdir,'GHCND_clean_complete_seasonal_tmin.txt'))
  rm(tmin_ghcn); gc()
  

# Get indices:
#------------------------------------------------------------------
source('Code/prep_station_data.R')

# for topowx:
tmax_topo_seas <- read.table('Data/Station_dataframes/Topowx_homog_complete_seasonal_tmax.txt',stringsAsFactors=F)
tic('get topowx tmax indices') # took about 50 minutes
out <- prep_station_data(tmax_topo_seas,save_dir=F)
toc()
write.table(out,paste0(outdir,'with_indices/Topowx_homog_tmax_master.txt'))
rm(tmax_topo_seas,out); gc()
        
tmin_topo_seas <- read.table('Data/Station_dataframes/Topowx_homog_complete_seasonal_tmin.txt',stringsAsFactors=F)
tic('get topowx tmin indices')
out <- prep_station_data(tmin_topo_seas)
toc()
write.table(out,paste0(outdir,'with_indices/Topowx_homog_tmin_master.txt'))
rm(tmin_topo_seas,out); gc()


# for ghcnd:
tmax_ghcn_seas <- read.table('Data/Station_dataframes/GHCND_clean_complete_seasonal_tmax.txt',stringsAsFactors=F)
names(tmax_ghcn_seas)[1] <- 'station_id'
tmax_ghcn_seas <- tmax_ghcn_seas %>% select(-stn_name)
tic('get ghcnd tmax indices')
out <- prep_station_data(tmax_ghcn_seas)
toc()
write.table(out,paste0(outdir,'with_indices/GHCND_clean_tmax_master.txt'))
rm(tmax_ghcn_seas,out); gc()

tmin_ghcn_seas <- read.table('Data/Station_dataframes/GHCND_clean_complete_seasonal_tmin.txt',stringsAsFactors=F)
names(tmin_ghcn_seas)[1] <- 'station_id'
tmin_ghcn_seas <- tmin_ghcn_seas %>% select(-stn_name)
tic('get ghcnd tmin indices')
out <- prep_station_data(tmin_ghcn_seas)
toc()
write.table(out,paste0(outdir,'with_indices/GHCND_clean_tmin_master.txt'))
rm(tmin_ghcn_seas,out); gc()                                    





  
  
  
  