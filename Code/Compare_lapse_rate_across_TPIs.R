# compare lapse rates calculated from different subsets fo stations based on topography
# first calculate lapse rates based on the subset of stations with TPI>0 and the subset with TPI<0 for each TPI radius
# then randomly sample subsets of the same size as above and calculate lapse rates.  Do this 1000 times for each TPI radius.
# calculate the mean lapse rate and standard deviation (across years) for each season for each sample
# compare the difference in mean and sd of lapse rates between the initial two subsets to the 
# differences in the randomly sampled subsets by calculating where in the cdf of randomly sampled differences
# the TPI difference falls.

buffer <- c('m200','m1000','m5000','m500')

pertab <- data.frame(matrix(ncol=4,nrow=(length(buffer)*4)))
names(pertab) <- c('tpi','season','dmean','dsd')
pertab$tpi <- rep(buffer,each=4)
pertab$season <- rep(c(sets$season),length(buffer))

for (bb in 1:length(buffer)){
  #thistpi<-quo(m200)
  thistpi <- buffer[bb]
  
  out <- season_means %>% filter(get(thistpi) >0) %>% group_by(season,year) %>%
    summarise(lapse_high=summary(lm(smtemp~elev))$coefficients[2]*1000)
  
  out2 <- season_means %>% filter(get(thistpi) <0) %>% group_by(season,year) %>%
    summarise(lapse_low=summary(lm(smtemp~elev))$coefficients[2]*1000)
  
  outs <- left_join(out,out2, by=c('season','year'))
  
  highlow <- outs %>% group_by(season) %>% summarise(mhigh=mean(lapse_high),mlow=mean(lapse_low),
                                          shigh=sd(lapse_high),slow=sd(lapse_low))
  da <- highlow %>% group_by(season) %>%
    summarise(dmean=mhigh-mlow, dsd=shigh-slow) %>% gather('metric','value',dmean:dsd)
  
  #outs <- outs %>% gather('lapse_type','lapse',lapse_high:lapse_low)
  #ggplot(outs) + geom_line(aes(x=year,y=lapse,col=lapse_type,group=lapse_type)) + facet_wrap(~season)
  
  
  diflapse <- data.frame(matrix(ncol=3,nrow=0))
  names(diflapse) <- c('season','dmean','dsd')
  xx <- which(names(tpi)==thistpi)
  for (ii in 1:1000){
    s1 <- unique(season_means$station_id) %>% sample(length(which(tpi[xx]>0)))
    set1 <- season_means %>% filter(station_id %in% s1) %>% group_by(season,year) %>%
      mutate(lapse=summary(lm(smtemp~elev))$coefficients[2]*1000) %>% group_by(season) %>%
      summarise(mlapse1=mean(lapse,na.rm=T),slapse1=sd(lapse,na.rm=T))
    set2 <- season_means %>% filter(!station_id %in% s1) %>% group_by(season,year) %>%
      mutate(lapse=summary(lm(smtemp~elev))$coefficients[2]*1000) %>% group_by(season) %>%
      summarise(mlapse2=mean(lapse,na.rm=T),slapse2=sd(lapse,na.rm=T))
    sets <- left_join(set1,set2,by='season')
    #sets %>% group_by(season) %>% summarise(dmean=mlapse1-mlapse2,dsd=slapse1-slapse2)
    
    diflapse <- bind_rows(diflapse, (sets %>% group_by(season) %>% summarise(dmean=mlapse1-mlapse2,dsd=slapse1-slapse2)))
  
  }
  
  dl <- diflapse %>% gather('metric','value',dmean:dsd)
  
  plot_names <- list('dmean'="Difference in mean value",'dsd'="Difference in standard deviation")
    plot_labeller <- function(variable,value){
      return(plot_names[value])
    }
  
  ggplot(dl, aes(season,(value))) + geom_boxplot() + facet_wrap(~metric, ncol=1, scales="free",labeller=plot_labeller) +
    geom_point(inherit.aes=F, data=da, aes(x=season,y=(value),col='red'),show.legend=F) +
    labs(y='C/km',x='Season',title=paste0('Comparing lapse rates across ',substr(thistpi,2,nchar(thistpi)),'m TPI classes'))
  #ggsave(paste0('Figures/dLapse_across_TPI_',thistpi,'.jpeg'))
  

  for(seas in 1:4){
    ii <- which(pertab$tpi==thistpi & pertab$season==unique(dl$season)[seas])
    percentile <- ecdf(dl$value[which(dl$metric=='dmean' & dl$season==unique(dl$season)[seas])])
    pertab$dmean[ii] <- percentile(da$value[which(da$metric=='dmean' & da$season==unique(dl$season)[seas])])
    
    percentile <- ecdf(dl$value[which(dl$metric=='dsd' & dl$season==unique(dl$season)[seas])])
    pertab$dsd[ii] <- percentile(da$value[which(da$metric=='dsd' & da$season==unique(dl$season)[seas])])
    
  }
} 
  
  #write.table(pertab, 'Figures/dLapse_percentiles.txt')
  
  # calculate the average percentile of the mean difference and sd difference (across all seasons, all tpis):
  pertab %>% summarise(mlapse=mean(dmean),slapse=mean(dsd))
  # calculate the mean percentiles across seasons:
  pertab %>% group_by(season) %>% summarise(mlapse = mean(dmean),slapse =mean(dsd))
  
  
  
  
