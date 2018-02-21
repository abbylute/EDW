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
  labs(y='', x='', title=paste0('TopoWx ',varname1,' stations (n=',nrow(tmax_meta),')\nwith at least ',min_pyrs*100,'% of years \nwith at least ',min_pdays*100,'% of data since ',substr(startdate,1,4)) )

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

# relate the seasonal lapse rate time series to covariates:

# compare near-surface seasonal lapse rates to free-air seasonal lapse rates:
library(ncdf4)
dirr <- 'Data/era-i-heights/'
files <- list.files(dirr)
year <- substr(files,8,11)

free_air <- data.frame()
for (yy in 1:length(files)){
  nc <- nc_open(paste0(dirr,files[yy]))
  lons <- ncvar_get(nc,'longitude')-360
  lats <- ncvar_get(nc,'latitude')
  times <- ncvar_get(nc,'time') #these are just monthly values (12 total)
  levels <- ncvar_get(nc,'level')
  
  d1 <- which(lons>sub_bbox[1] & lons<sub_bbox[3])
  d2 <- which(lats>sub_bbox[2] & lats<sub_bbox[4])
  temp <- ncvar_get(nc,'t',start=c(min(d1),min(d2),1,1),count=c(length(d1),length(d2),-1,-1)) -273.15 # grab temperatures and convert to C
  geop <- ncvar_get(nc,'z',start=c(min(d1),min(d2),1,1),count=c(length(d1),length(d2),-1,-1))
  
  temp <- apply(temp,c(3,4),mean) # average over the bounding box to get a temperature for each height and month in this year
  geop <- apply(geop,c(3,4),mean)/9.8
  
  temp <- cbind(levels,data.frame(temp)); #names(temp2) <- c('height',paste0('m',c(1:12)))
  temp <- gather(temp,key=month,value=temperature, 2:13)
  temp$month <- factor(temp$month)
  levels(temp$month)=c('01','02','03','04','05','06','07','08','09','10','11','12')
  temp <- temp %>% left_join(seasons, by = "month")
  
  geop <- cbind(levels,data.frame(geop)); #names(temp2) <- c('height',paste0('m',c(1:12)))
  geop <- gather(geop,key=month,value=height, 2:13)
  geop$month <- factor(geop$month)#, 
  levels(geop$month)=c('01','02','03','04','05','06','07','08','09','10','11','12')
  geop <- geop %>% left_join(seasons, by='month')
  
  jj <- left_join(temp,geop, by=c('levels','month','season'))
  jj <- jj %>% group_by(season) %>% summarise(lapse=summary(lm(temperature~height))$coefficients[2]*1000)
  jj <- cbind(year=rep(year[yy],4),jj, stringsAsFactors=F)
  
  free_air <- rbind(free_air, jj)
} # end yearly file loop



unique(free_air$year)
unique(season_lapse_rates$year)

fa <- free_air %>% filter(year>1979, year<2016) 
ns <- season_lapse_rates %>% filter(year>1979, year<2016)
fa <- select(fa,year,season,lapse); names(fa) <- c('year','season','fa_lapse')
ns <- select(ns,year,season,lapse); names(ns) <- c('year','season','ns_lapse')
lapse <- join(fa,ns)

lapse %>% group_by(season) %>%
  summarise(corr=cor(fa_lapse,ns_lapse))

lapse2 <- lapse %>% gather(key='metric',value='lapse',fa_lapse,ns_lapse)

ggplot(lapse2,aes(x=year,y=lapse,group=metric,color=metric)) +
  geom_line() +
  facet_wrap(~season)






