library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)

varname <- 'tmax'  # specify tmax or tmin
meta <- tmax_meta
temps <- tmax

sub_bbox <- c(-113,42,-110,45)

ss <- which(meta$lat>sub_bbox[2] & meta$lat<sub_bbox[4] & meta$lon>sub_bbox[1] & meta$lon<sub_bbox[3])
meta<- meta[ss,]
temps <- data.frame(temps[ss,])  
cal <- seq(as.Date(startdate),as.Date(enddate),'days')

#ggmap(gm) +
#  geom_point(data=sub_tmax_meta, aes(x=lon, y=lat,color=elev),size=1)+
##  scale_color_gradientn(name="Elevation", colors=rainbow(20))+
#  theme(legend.key=element_rect(fill="white")) +
#  coord_map(ylim=bbox[c(2,4)], xlim=bbox[c(1,3)]) +
#  labs(y='', x='' )

temps <- cbind(meta,temps)
tm <- temps %>% gather(key=date, value=temperature, 6:ncol(temps))

datedf <- data.frame(date=c(colnames(temps[6:ncol(temps)])), datenm=cal, stringsAsFactors=F)
tm <- left_join(tm,datedf, by='date')
tm <- tm %>% select(-date)
names(tm)[ncol(tm)] <- 'date'
tm <- tm %>% separate(date, into=c('year','month','day'),sep='-',remove=F)
tm <- arrange(tm,station_id,date)


seasons <- data.frame('Fall'=c('09','10','11'), 'Winter'=c('12','01','02'),'Spring'=c('03','04','05'),'Summer'=c('06','07','08'))
  seasons <- gather(seasons,season,month)

tm2 <- tm %>% left_join(seasons, by = "month")

season_means <- tm2 %>%
  group_by(station_id,season,year) %>%
  mutate(smtemp = mean(temperature,na.rm=T)) %>% select(-temperature,-date,-month,-day)
season_means <- season_means %>% group_by(station_id,season,year) %>% slice(which.min(smtemp))

season_lapse_rates <- season_means %>% group_by(season,year) %>%
  summarise(lapse=summary(lm(smtemp~elev))$coefficients[2]*1000)

ggplot(data=season_lapse_rates,aes(x=year,y=lapse,group=season)) +  geom_line() + facet_wrap(~season)


#################
# Compute free-air lapse rates for each lat/lon column for each year and month using only levels above the reanalysis land surface
rtse= T  # if rtse (restrict to station elevations) is T, then only free-air observations within the range of stations elevations are used to calculate free-air lapse rates
         # if rtse is F, then free-air observations are only bounded on the bottom by the reanalysis land surface height
         # NOTE: this removes a lot of observations since many free air obs are at >3000m
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
  
  if (rtse==T){ # to restrict free-air lapse rate calculations to the range of station elevations for comparison
    dd <- dd %>% filter(geop > (min(meta$elev)-200) & geop < (max(meta$elev)+200))
  }
  
  ## Calculate lapse rate in each lat/lon column
  dd <- dd %>% group_by(lon,lat,month) %>%
    mutate(lapse = summary(lm(temp~geop))$coefficients[2]*1000)

  dd <- bind_cols(year=rep(year[yy],nrow(dd)),dd)
  free_air <- bind_rows(free_air, dd)
} # end yearly file loop

# join free_air df with seasons to make it compatible with the station lapse rate calculations:
  seasons_numeric <- seasons; seasons_numeric$month <- as.integer(seasons_numeric$month)
  free_air <- join(free_air,seasons_numeric, by='month')
  
# calculate the seasonal mean lapse rates as the average of the monthly lapse rates already calculated
  fa <- free_air %>% group_by(year,lat,lon,season) %>% summarise(lapse_seasonal = mean(lapse, na.rm=T))
# average seasonal lapse rates over the domain for comparison with stations:
  fa <- fa %>% group_by(year,season) %>% summarise(lapse_seasonal=mean(lapse_seasonal, na.rm=T))

  
# join the near-surface and free-air lapse rates:
  unique(free_air$year)
  unique(season_lapse_rates$year)
  fa <- fa %>% filter(year>1979, year<2016) 
  ns <- season_lapse_rates %>% filter(year>1979, year<2016)
  names(fa) <- c('year','season','fa_lapse')
  ns <- select(ns,year,season,lapse); names(ns) <- c('year','season','ns_lapse')
  lapse <- left_join(fa,ns, by=c('year','season'))

  
# seasonal correlations between free-air and near-surface lapse rates:
corr <- lapse %>% group_by(season) %>%
  summarise(corr=signif(cor(fa_lapse,ns_lapse)^2,2)) # this gives R-squared value


lapse2 <- lapse %>% gather(key='metric',value='lapse',fa_lapse,ns_lapse)

gg <- ggplot(lapse2,aes(x=year,y=lapse,group=metric,color=metric)) +
  geom_line() +
  geom_text(inherit.aes=F, data=corr, aes(x=6,y=-1,label=paste0('R^2=',corr))) +
  labs(y='Lapse Rate (C/km)', title=paste0('Seasonal ',varname,' lapse rates'),subtitle=substr(figdir,9,nchar(figdir)-1)) +
  scale_color_manual('',values=c('darkgreen','purple'),labels=c('free-air','near-surface')) +
  scale_x_discrete(labels=seq(1980,2020,5),breaks=seq(1980,2020,5)) +
  facet_wrap(~factor(season,levels=c("Spring","Summer","Fall","Winter")))
gg

jpeg(filename=paste0(figdir,varname,'_seasonal_lapse_rates.jpeg'),width=10,height=7,units="in",res=500,quality=100)
print(gg)
dev.off()


#gg <- ggplot(lapse,aes(x=ns_lapse,y=fa_lapse)) +
#  geom_point() +
#  facet_wrap(~season)
#gg
#xx <- which(lapse$season=='Summer')
##plot(lapse$ns_lapse[xx],lapse$fa_lapse[xx])
#abline(0,1)

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

