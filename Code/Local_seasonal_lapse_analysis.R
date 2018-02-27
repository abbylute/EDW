library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)

sub_bbox <- c(-113,42,-110,45)

ss <- which(tmax_meta$lat>sub_bbox[2] & tmax_meta$lat<sub_bbox[4] & tmax_meta$lon>sub_bbox[1] & tmax_meta$lon<sub_bbox[3])
sub_tmax_meta<- tmax_meta[ss,]
sub_tmax <- data.frame(tmax[ss,])  
cal <- seq(as.Date(startdate),as.Date(enddate),'days')

ggmap(gm) +
  geom_point(data=sub_tmax_meta, aes(x=lon, y=lat,color=elev),size=1)+
  scale_color_gradientn(name="Elevation", colors=rainbow(20))+
  theme(legend.key=element_rect(fill="white")) +
  coord_map(ylim=bbox[c(2,4)], xlim=bbox[c(1,3)]) +
  labs(y='', x='' )

sub_tmax <- cbind(sub_tmax_meta,sub_tmax)
tm <- sub_tmax %>% gather(key=date, value=tmax, 6:ncol(sub_tmax))
tm$date <- factor(tm$date)
levels(tm$date) <- cal
tm <- tm %>% separate(date, into=c('year','month','day'),sep='-',remove=F)
tm <- arrange(tm,station_id,date)


seasons <- data.frame('Fall'=c('09','10','11'), 'Winter'=c('12','01','02'),'Spring'=c('03','04','05'),'Summer'=c('06','07','08'))
  seasons <- gather(seasons,season,month)

tm2 <- tm %>% left_join(seasons, by = "month")

season_means <- tm2 %>%
  group_by(station_id,season,year) %>%
  mutate(smtmax = mean(tmax,na.rm=T)) %>% select(-tmax,-date,-month,-day)
season_means <- season_means %>% group_by(station_id,season,year) %>% slice(which.min(smtmax))


season_lapse_rates <- season_means %>% group_by(season,year) %>%
  summarise(lapse=summary(lm(smtmax~elev))$coefficients[2]*1000)

ggplot(data=season_lapse_rates,aes(x=year,y=lapse,group=season)) +  geom_line() + facet_wrap(~season)


#################
# Compute free-air lapse rates for each lat/lon column for each year and month using only levels above the reanalysis land surface
library(ncdf4)
  dirr <- 'Data/era-i-heights/'
  files <- list.files(dirr, pattern='heights')
  year <- substr(files,8,11)
  
  # get the land surface heights to restrict the height range over which lapse rates are calculated:
  land <- nc_open(paste0(dirr,'era_land.nc'),write=F)  # contains a 'z' value for land surface height for each lat/lon. need to divide this by 9.8 as well.
  lons <- ncvar_get(land,'longitude')-360
  lats <- ncvar_get(land,'latitude')
  d1 <- which(lons>sub_bbox[1] & lons<sub_bbox[3])
  d2 <- which(lats>sub_bbox[2] & lats<sub_bbox[4])
  latlookup <- data.frame(lat =c(paste0('lat',c(1:length(d2)))), latval=lats[d2], stringsAsFactors=F)
  landheight <- ncvar_get(land,'z',start=c(min(d1),min(d2),1),count=c(length(d1),length(d2),-1))/9.8 # now in meters
  lh <- data.frame(lons[d1],landheight); names(lh) <- c('lon','lat1','lat2','lat3')
  lh <- gather(lh,key=lat, value=landz, lat1:lat3)
  lh <- lh %>% left_join(latlookup, by="lat") %>% select(lon,latval,landz); names(lh) <- c('lon','lat','landz')
  
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
  #########
  # make temp into a data.frame:
  tt <- as.data.frame.table(temp,responseName='temp')
  names(tt) <- c('lon','lat','level','month','temp')
  
  # make geop into a data.frame:
  dd <- as.data.frame.table(geop,responseName='geop')
  names(dd) <- c('lon','lat','level','month','geop')
  dd <- join(dd,tt, by=c('lon','lat','level','month')) # join in the temperature data
  dd <- dd %>% left_join(lonlookup, by="lon") %>% left_join(latlookup, by="lat") %>% 
    left_join(levlookup, by="level") %>% left_join(timelookup, by="month") 
  dd <- dd %>% select(lonval:timeval,geop,temp);   names(dd) <- c('lon','lat','level','month','geop','temp')

  dd <- join(dd,lh, by=c('lon','lat')) # join with land height data
  dd <- dd %>% filter(geop > landz) # remove observations where the geopotential height measurement is below land surface
  
  ##########
  #temp <- apply(temp,c(3,4),mean) # average over the bounding box to get a temperature for each height and month in this year
  #geop <- apply(geop,c(3,4),mean)
  
  ## Calculate lapse rate in each lat/lon column
  dd <- dd %>% group_by(lon,lat,month) %>%
    mutate(lapse = summary(lm(temp~geop))$coefficients[2]*1000)
  
  
  #temp <- cbind(levels,data.frame(temp)); #names(temp2) <- c('height',paste0('m',c(1:12)))
  #temp <- gather(temp,key=month,value=temperature, 2:13)
  #temp$month <- factor(temp$month)
  #levels(temp$month)=c('01','02','03','04','05','06','07','08','09','10','11','12')
  #temp <- temp %>% left_join(seasons, by = "month")
  
  #geop <- cbind(levels,data.frame(geop)); #names(temp2) <- c('height',paste0('m',c(1:12)))
  #geop <- gather(geop,key=month,value=height, 2:13)
  #geop$month <- factor(geop$month)#, 
  #levels(geop$month)=c('01','02','03','04','05','06','07','08','09','10','11','12')
  #geop <- geop %>% left_join(seasons, by='month')
  
  #jj <- left_join(temp,geop, by=c('levels','month','season'))
  #jj <- jj %>% group_by(season) %>% summarise(lapse=summary(lm(temperature~height))$coefficients[2]*1000)
  #jj <- jj %>% group_by(month,season) %>% summarise(lapse=summary(lm(temperature~height))$coefficients[2]*1000)
  dd <- bind_cols(year=rep(year[yy],nrow(dd)),dd)
  #jj <- cbind(year=rep(year[yy],4),jj, stringsAsFactors=F)
  free_air <- bind_rows(free_air, dd)
} # end yearly file loop



unique(free_air$year)
unique(season_lapse_rates$year)

fa <- free_air %>% filter(year>1979, year<2016) 
ns <- season_lapse_rates %>% filter(year>1979, year<2016)
fa <- select(fa,year,season,lapse); names(fa) <- c('year','season','fa_lapse')
ns <- select(ns,year,season,lapse); names(ns) <- c('year','season','ns_lapse')
lapse <- join(fa,ns)

lapse %>% group_by(season) %>%
  summarise(corr=cor(fa_lapse,ns_lapse)^2) # this gives R-squared value

lapse2 <- lapse %>% gather(key='metric',value='lapse',fa_lapse,ns_lapse)

ggplot(lapse2,aes(x=year,y=lapse,group=metric,color=metric)) +
  geom_line() +
  facet_wrap(~season)




#################
#TRYING TO OPEN CALIPSO DATA:
library(rgdal)
library(gdalUtils)
fn <- 'Documents/IGERT/EDW/DATA/CALIPSO/CAL_LID_L2_05kmALay-Standard-V4-10.2006-06-24T20-12-50ZD_Subset.hdf'
out <- readGDAL(fn, drivername='HDF5')
gdalinfo(fn)
subfn <- get_subdatasets(fn)
gdal_translate(subfn[1],'Documents/IGERT/EDW/DATA/CALIPSO/output.tif')

library(ncdf4)
nc <- nc_open(fn)

library(raster)
r <- raster('Documents/IGERT/EDW/DATA/CALIPSO/output1.tif')
r
hist(r)
image(r)


aa <- data.frame(rasterToPoints(r))

aa <- data.frame(rasterToPoints(aggregate(r, fact=c(1,dim(r)[1]) , fun=fun, na.rm=T)))

