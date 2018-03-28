# calculate free air lapse rates from reanalysis data:
# bbox is the bounding box within which to compute free-air lapse rates
# rtse (restrict to station elevations): if rtse is T (default), then only free-air observations within the range of stations elevations are used to calculate free-air lapse rates
# if rtse is F, then free-air observations are only bounded on the bottom by the reanalysis land surface height
# NOTE: this removes a lot of observations since many free air obs are at >3000m


get_free_air_lapse_rates = function(bbox, rtse=T){
  
  library(ncdf4)
  dirr <- 'Data/era-i-heights/'
  files <- list.files(dirr, pattern='heights')
  year <- substr(files,8,11)
  
  # get the land surface heights to restrict the height range over which lapse rates are calculated:
  land <- nc_open(paste0(dirr,'era_land.nc'),write=F)  # contains a 'z' value for land surface height for each lat/lon. need to divide this by 9.8 as well.
  lons <- ncvar_get(land,'longitude')-360
  lats <- ncvar_get(land,'latitude')
  d1 <- which(lons>bbox[1] & lons<bbox[3])
  d2 <- which(lats>bbox[2] & lats<bbox[4])
  latlookup <- data.frame(lat =c(paste0('lat',c(1:length(d2)))), latval=lats[d2], stringsAsFactors=F)
  landheight <- ncvar_get(land,'z',start=c(min(d1),min(d2),1),count=c(length(d1),length(d2),-1))/9.8 # now in meters
  lh <- data.frame(lons[d1],landheight); names(lh) <- c('lon','lat1','lat2','lat3')
  lh <- gather(lh,key=lat, value=landz, lat1:lat3)
  lh <- lh %>% left_join(latlookup, by="lat") %>% dplyr::select(lon,latval,landz); names(lh) <- c('lon','lat','landz')

  free_air <- data.frame()
  for (yy in 1:length(files)){
    nc <- nc_open(paste0(dirr,files[yy]))
    lons <- ncvar_get(nc,'longitude')-360
    lats <- ncvar_get(nc,'latitude')
    times <- ncvar_get(nc,'time') #these are just monthly values (12 total)
    levels <- ncvar_get(nc,'level')
    
    d1 <- which(lons>sub_bbox[1] & lons<sub_bbox[3])
    d2 <- which(lats>sub_bbox[2] & lats<sub_bbox[4])
    latlookup <- data.frame(lat =c(paste0(LETTERS[1:length(d2)])), latval=lats[d2])
    lonlookup <- data.frame(lon =c(paste0(LETTERS[1:length(d1)])), lonval=lons[d1])
    timelookup <- data.frame(month=c(paste0(LETTERS[1:length(times)])), timeval=c(1:12))
    levlookup <- data.frame(level =c(paste0(LETTERS[1:length(levels)])), levval=levels)
    
    temp <- ncvar_get(nc,'t',start=c(min(d1),min(d2),1,1),count=c(length(d1),length(d2),-1,-1)) -273.15 # grab temperatures and convert to C
    geop <- ncvar_get(nc,'z',start=c(min(d1),min(d2),1,1),count=c(length(d1),length(d2),-1,-1)) / 9.8
    
    # make temp into a data.frame:
    tt <- as.data.frame.table(temp,responseName='temp')
    names(tt) <- c('lon','lat','level','month','temp')
    
    # make geop into a data.frame:
    dd <- as.data.frame.table(geop,responseName='geop')
    names(dd) <- c('lon','lat','level','month','geop')
    dd <- join(dd,tt, by=c('lon','lat','level','month')) # join in the temperature data
    dd <- dd %>% left_join(lonlookup, by="lon") %>% left_join(latlookup, by="lat") %>% 
      left_join(levlookup, by="level") %>% left_join(timelookup, by="month") 
    dd <- dd %>% dplyr::select(lonval:timeval,geop,temp);   names(dd) <- c('lon','lat','level','month','geop','temp')
    
    dd <- join(dd,lh, by=c('lon','lat')) # join with land height data
    dd <- dd %>% filter(geop > landz) # remove observations where the geopotential height measurement is below land surface
    
    if (rtse==T){ # to restrict free-air lapse rate calculations to the range of station elevations for comparison
      dd <- dd %>% filter(geop > (min(meta$elev)-200) & geop < (max(meta$elev)+200))
    }
    
    ## Calculate lapse rate in each lat/lon column
    dd <- dd %>% group_by(lon,lat,month) %>%
      mutate(lapse = summary(lm(temp~geop))$coefficients[2]*1000)
    
    dd <- bind_cols(year=rep(year[yy],nrow(dd)),dd)
    free_air <- bind_rows(free_air, dd)
  } # end yearly file loop
  return(free_air)
}