library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggmap)
source('Code/FIX_get_stamenmap_RCode.R') #ignore png warning

#-----------------------------------------------------------------      
# SETUP:
master <- read.table("Data/TopoWx/GHCN_raw_1980_summary.txt")
regname <- 'Upper Snake'
figdir <- paste0("Figures/GHCN_raw_1980/",gsub(' ','',regname),"/")
sub_bbox <- c(-113,42,-110,45)
#sub_bbox <- c(-114,41,-109,46) # wider bbox
#sub_bbox <- c(-122,41,-118,46) # east side of oregon cascades
#sub_bbox <- c(-110,37,-105,41) # western colorado


df <- master %>% filter(lat>sub_bbox[2] & lat<sub_bbox[4] & lon>sub_bbox[1] & lon<sub_bbox[3])
length(unique(df$station_id))

cal <- seq(as.Date(startdate),as.Date(enddate),'days')

#-----------------------------------------------------------------      
# Map the selected stations:
gm <- get_stamenmap(sub_bbox,maptype="toner",zoom=6)  
ggmap(gm) +
  geom_point(data=df, aes(x=lon, y=lat,color=elev),size=1)+
  scale_color_gradientn(name="Elevation", colors=rainbow(20))+
  theme(legend.key=element_rect(fill="white")) +
  coord_map(ylim=sub_bbox[c(2,4)], xlim=sub_bbox[c(1,3)]) +
  labs(y='', x='' )

#-----------------------------------------------------------------      
## VARIABLE CORRELATIONS:
tvars <- c('tmax','tmin')
exvars <- c('elev','srad','m200','m500','m1000','m5000')
seas <- unique(df$season)

tabm <- df %>% group_by(season,station_id) %>% 
  summarise(elev=mean(elev,na.rm=T),tmax=mean(tmax,na.rm=T),tmin=mean(tmin,na.rm=T),
            srad=mean(srad,na.rm=T),m200=mean(m200,na.rm=T),m500=mean(m500,na.rm=T),m1000=mean(m1000,na.rm=T),m5000=mean(m5000,na.rm=T))
cortab <- data.frame('tvar' = rep(tvars, each=length(exvars)*length(seas)),
                     'var' = rep(rep(exvars, each=length(seas)), length(tvars)),
                     'season' = rep(seas, length(exvars)*length(tvars)),
                     'cor' = rep(NA, length(tvars)*length(exvars)*length(seas)))

for (tt in 1:length(tvars)){
  for (ss in 1:length(seas)){
    tabb <- tabm %>% filter(season==seas[ss]) %>% select(tvars[tt],station_id,season,exvars) %>% na.omit
    for (vv in 1:length(exvars)){
      rr <- which(cortab$tvar==tvars[tt] & cortab$var==exvars[vv] & cortab$season==seas[ss])
      cortab$cor[rr] <- cor(tabb[which(names(tabb)==tvars[tt])], tabb[which(names(tabb)==exvars[vv])])
    }
  }
}
cortab$var <- factor(cortab$var, unique(cortab$var))
cortab %>% group_by(tvar,var) %>% summarise(mcor=mean(cor))

ggplot(cortab) + geom_line(aes(x=season,y=cor,col=var,group=var)) + facet_wrap(~tvar, nrow=3) +
  scale_color_discrete('Explanatory\nVariable', labels=c('elevation','solar radiation','TPI 200m','TPI 500m','TPI 1km','TPI 5km')) +
  labs(y="correlation (R, across sites)", title=paste0(regname,' TopoWx temperature-variable correlations'))
ggsave(paste0(figdir,'temperature_variable_correlations.jpeg'))


# LAPSE RATE METHOD 1: Calculate seasonal lapse rates using elevation alone:
#-------------------------------------------------------------------------- 
  season_means <- season_means %>% group_by(season,year) %>%
  mutate(lapse_elev=summary(lm(smtemp~elev))$coefficients[2]*1000)
  #ggplot(data=season_means,aes(x=year,y=lapse_elev,group=season)) +  geom_line() + facet_wrap(~season)
  #season_means %>% group_by(season) %>% summarise(ml=mean(lapse_elev))
  

# LAPSE RATE METHOD 2: Calculate seasonal lapse rates using elevation and TPI:
#-------------------------------------------------------------------------- 
  buffer <- names(tpi)[2:ncol(tpi)] # need to update this so that buffer comes from teh column names
    
  for (ii in 1:length(buffer)){
    var_name <-paste0('lapse_',buffer[ii])
      season_means <- season_means %>% group_by(season,year) %>%
                       mutate(!!var_name :=summary(lm(as.formula(paste0('smtemp ~ elev + ', buffer[ii]))))$coefficients[2]*1000)
  }

  #season_means_long <- season_means %>% 
  #  gather('lapse_type','lapses',c('lapse_elev',paste0('lapse_',buffer)))
  #season_means_long$lapse_type <- factor(season_means_long$lapse_type,unique(season_means_long$lapse_type))
  #ggplot(data=season_means_long,aes(x=year,y=lapses,group=lapse_type,col=lapse_type)) +  geom_line() + facet_wrap(~season)
  
  
# LAPSE RATE METHOD 3: Calculate seasonal lapse rates using elevation and solar radiation:
#-------------------------------------------------------------------------- 
  # calculate lapse rates based on elevation and seasonal solar radiation:
  season_means <- season_means %>% group_by(season,year) %>%
    mutate(lapse_srad=summary(lm(smtemp~elev+srad))$coefficients[2]*1000)
  
  season_means_long <- season_means %>% 
    gather('lapse_type','lapses',c('lapse_elev',paste0('lapse_',buffer),'lapse_srad'))
  season_means_long$lapse_type <- factor(season_means_long$lapse_type,unique(season_means_long$lapse_type))
  ggplot(data=season_means_long,aes(x=year,y=lapses,group=lapse_type,col=lapse_type)) +  geom_line() + facet_wrap(~season)
  
  # incorporation of solar radiation results in even more extreme winter lapse rates!
  ggplot(data=season_means) +
    geom_point(aes(x=lapse_elev,y=lapse_srad,col=season)) +
    geom_abline(aes(slope=1,intercept=0)) +
    labs(x='Lapse rate (elevation only)', y='Lapse rate (elevation and solar radiation')
  ggsave(paste0(figdir,'solar_vs_elev_lapse_rate.jpeg'))
  
  
# LAPSE RATE METHOD 4: Calculate seasonal lapse rates using elevation, TPI, and solar radiation:
#-------------------------------------------------------------------------- 
  for (ii in 1:length(buffer)){
    var_name <-paste0('lapse_srad_',buffer[ii])
    season_means <- season_means %>% group_by(season,year) %>%
      mutate(!!var_name :=summary(lm(as.formula(paste0('smtemp ~ elev + srad +', buffer[ii]))))$coefficients[2]*1000)
  }
  #season_means_long <- season_means %>% 
  #  gather('lapse_type','lapses',c('lapse_elev',paste0('lapse_',buffer),'lapse_srad',paste0('lapse_srad_',buffer)))
  #season_means_long$lapse_type <- factor(season_means_long$lapse_type,unique(season_means_long$lapse_type))
  #ggplot(data=season_means_long,aes(x=year,y=lapses,group=lapse_type,col=lapse_type)) +  geom_line() + facet_wrap(~season)
  
  

# Compute FREE-AIR LAPSE RATES for each lat/lon column for each year and month using only levels above the reanalysis land surface
#-------------------------------------------------------------------------- 
  source('Code/get_free_air_lapse_rates.R')
  free_air <- get_free_air_lapse_rates(sub_bbox,rtse=T,meta)
  
  # join free_air df with seasons to make it compatible with the station lapse rate calculations:
    seasons_numeric <- seasons; seasons_numeric$month <- as.integer(seasons_numeric$month)
    free_air <- join(free_air,seasons_numeric, by='month')
  
  # calculate the seasonal mean lapse rates as the average of the monthly lapse rates already calculated
    free_air <- free_air %>% group_by(year,lat,lon,season) %>% summarise(lapse_seasonal = mean(lapse, na.rm=T))
  # average seasonal lapse rates over the domain for comparison with stations:
    free_air <- free_air %>% group_by(year,season) %>% summarise(lapse=mean(lapse_seasonal, na.rm=T))

  # join the near-surface and free-air lapse rates:
    free_air <- free_air %>% filter(year>=min(season_means$year), year<=max(season_means$year)) 
    names(free_air) <- c('year','season','fa_lapse')
    season_means <- left_join(season_means, free_air,by=c('year','season'))

    season_means_long <- season_means %>% 
      gather('lapse_type','lapses',c('lapse_elev',paste0('lapse_',buffer),'lapse_srad',paste0('lapse_srad_',buffer),'fa_lapse'))
    season_means_long$lapse_type <- factor(season_means_long$lapse_type,unique(season_means_long$lapse_type))
    ggplot(data=season_means_long,aes(x=year,y=lapses,group=lapse_type,col=lapse_type)) +  geom_line() + facet_wrap(~season)
    
    
###################  
# Calculate seasonal correlations between free-air and all near-surface lapse rates:
fa_cors <- season_means %>% group_by(season) %>%
  summarise(fa_elev=signif(cor(fa_lapse,lapse_elev)^2,2),
            fa_m200=signif(cor(fa_lapse,lapse_m200)^2,2),
            fa_m5000=signif(cor(fa_lapse,lapse_m5000)^2,2),
            fa_m1000=signif(cor(fa_lapse,lapse_m1000)^2,2),
            fa_m50000=signif(cor(fa_lapse,lapse_m50000)^2,2),
            fa_srad=signif(cor(fa_lapse,lapse_srad)^2,2),
            fa_srad_m200=signif(cor(fa_lapse,lapse_srad_m200)^2,2),
            fa_srad_m5000=signif(cor(fa_lapse,lapse_srad_m5000)^2,2),
            fa_srad_m1000=signif(cor(fa_lapse,lapse_srad_m1000)^2,2),
            fa_srad_m50000=signif(cor(fa_lapse,lapse_srad_m50000)^2,2)
            ) # this gives R-squared value
    write.table(fa_cors,paste0(figdir,'corr_freeair_nslapse_methods.txt'))
  cbind(fa_cors$season,rowMeans(fa_cors[,2:ncol(fa_cors)])) 
  # in general, there are stronger correlations between free-air and near-surface lapse methods during fall and winter
  colMeans(fa_cors[,2:ncol(fa_cors)])
  # averaged across seasons, the near-surface lapse method most strongly correlated with free-air is the one just using elevation!  
    
  ntype <- length(unique(season_means_long$lapse_type))
  gg <- ggplot(season_means_long,aes(x=year,y=lapses,group=lapse_type,color=lapse_type)) +
    geom_line() +
    #geom_text(inherit.aes=F, data=corr, aes(x=6,y=-1,label=paste0('R^2=',corr))) +
    labs(y='Lapse Rate (C/km)', title=paste0('Seasonal ',varname,' lapse rates'),subtitle=substr(figdir,9,nchar(figdir)-1)) +
    scale_color_manual('',values=c(rainbow(ntype)[-ntype],'black'),  #c('darkgreen','purple'),
                       labels=c(gsub('lapse_','',unique(season_means_long$lapse_type))[-ntype],'free-air')) +
    scale_x_discrete(labels=seq(1980,2020,5),breaks=seq(1980,2020,5)) +
    facet_wrap(~factor(season,levels=c("Spring","Summer","Fall","Winter")))
  gg

  jpeg(filename=paste0(figdir,varname,'_seasonal_lapse_rates.jpeg'),width=10,height=7,units="in",res=500,quality=100)
  print(gg)
  dev.off()




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

