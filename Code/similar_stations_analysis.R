# Similar Stations Analysis:


similar_stations_analysis = function(data, vars, opt, tvar, figdir){
  
  if (!all(c('station_id',vars,tvar,'year','season') %in% names(data))){
    print('Error: dataframe called meta must contain column names station_id and the names of the specified variables to use for calculating similarity (vars)')
  }
  varnames <- paste0('var',c(1:length(vars)))
  names(data)[which(names(data)==tvar)] <- 'tvar'
  if ('elev' %in% vars){
    data <- cbind(data,'elev2'=data$elev)
  }
  names(data)[which(names(data) %in% vars)] <- varnames
  if ('elev2' %in% names(data)){
    names(data)[which(names(data) == 'elev2')] <- 'elev'
  }
  
  
#  sm <- data %>% dplyr::select(station_id,tvar,varnames) %>% 
 #   group_by(station_id) %>% summarise_(tvar=mean(tvar),  varnames=mean(eval(parse(text = varnames))))
  sm <- data %>% group_by(station_id, add=F) %>% dplyr::select(station_id,tvar,varnames)
  sm <- sm %>% group_by(station_id) %>% summarise_all(mean)
  
  #sm <- data %>% dplyr::select(station_id,tvar,varnames) %>% 
  #  group_by(station_id) %>% summarise(tvar=mean(tvar),  eval(parse(text = varnames))=mean(eval(parse(text = varnames))))
  #,msrad=mean(srad)) %>% 
  #  slice(which.min(mtemp)) %>% dplyr::select(-c(year,vars,season,tvar))
  
  #sm <- data %>% dplyr::select(station_id,tvar, vars) %>% 
  #  group_by(station_id) %>% mutate(mtemp=mean(smtemp),msrad=mean(srad)) %>% 
  #  slice(which.min(mtemp)) %>% dplyr::select(-c(year,smtemp,season,srad))
  # note that whether you scale before aggregating across years/seasons or after 
  # makes a big difference in the scaled values, and the stations picked.
  
  source('Code/get_similar_stations.R')
  
  samp_size <- c(2:nrow(sm))
  lapse_tab <- data.frame('year'=rep(unique(data$year),length(unique(data$season))), 'season'=rep(unique(data$season),each=length(unique(data$year))))
  lapse_tab_rmse <- lapse_tab
  picks <- list()
  elevrng <- numeric()
  for (nn in samp_size){ # for this range of sample sizes, should include the total sample size
    picks[[nn-1]] <- get_similar_stations(sm, vars=varnames, opt=opt,num_stations=nn, scale=T)
    
    lapse_tab <- left_join(lapse_tab, data %>% filter(station_id %in% picks[[nn-1]]) %>%
                             group_by(year,season, add=F) %>%  # add=F removes the site grouping
                             summarise(samp_lapse = summary(lm(tvar~elev))$coefficients[2]*1000),
                           by=c('year','season'))
    lapse_tab_rmse <- left_join(lapse_tab_rmse, data %>% filter(station_id %in% picks[[nn-1]]) %>%
                                  group_by(year,season, add=F) %>%  # add=F removes the site grouping
                                  summarise(lapse_rmse = sqrt(sum(summary(lm(tvar~elev))$residuals^2)/nn)),
                                by=c('year','season'))
    names(lapse_tab)[ncol(lapse_tab)] <- paste0('nn',nn)
    names(lapse_tab_rmse)[ncol(lapse_tab_rmse)] <- paste0('nn',nn)
    
    elevrng[nn-1] <- as.numeric(data %>% filter(station_id %in% picks[[nn-1]]) %>% group_by() %>% summarise(elevrng =diff(range(elev))))
  }
  lapse_tab <- lapse_tab %>% gather('num','lapse',paste0('nn',min(samp_size)):paste0('nn',max(samp_size)))
  lapse_tab <- left_join(lapse_tab, data.frame(cbind(elevrng, 'num'=paste0('nn',2:nrow(sm)))), by='num')
  lapse_tab$num <- factor(lapse_tab$num, unique(lapse_tab$num)) # this is very important to display the results correctly!
  
  lapse_tab_rmse <- lapse_tab_rmse %>% gather('num','lapse_rmse',paste0('nn',min(samp_size)):paste0('nn',max(samp_size)))
  lapse_tab_rmse <- left_join(lapse_tab_rmse, data.frame(cbind(elevrng, 'num'=paste0('nn',2:nrow(sm)))), by='num')
  lapse_tab_rmse$num <- factor(lapse_tab_rmse$num, unique(lapse_tab_rmse$num)) # this is very important to display the results correctly!
  ## NOTE ## Sometimes a station that is selected may be missing data for one or more years
  
  lapse_tab <- left_join(lapse_tab, lapse_tab_rmse, by=c('year','season','num','elevrng'))
  
  
  ggplot(lapse_tab) +geom_line(aes(x=year,y=lapse,group=num,col=num)) + 
    labs(y='Lapse Rate (C/km)', title=paste0(tvar,' lapse rates for samples based on ', paste(vars,collapse=', '))) +
    scale_color_manual(name='sample size',labels=seq(2,nrow(sm),1),values=c(rainbow(length(samp_size)-1),'black')) +
    facet_wrap(~season)# +
    #coord_cartesian(ylim=c(-20,20))
  #coord_cartesian(ylim=c(pmax(-14,min(lapse_tab$lapse)), pmin(10,max(lapse_tab$lapse))))
  ggsave(paste0(figdir,tvar,'_similar_station_seasonal_lapse_rates.jpeg'))
  
  lapse_tab %>% group_by(season,num) %>% summarise(sdl = sd(lapse)) %>%
    ggplot() + geom_point(aes(x=num,y=sdl)) +facet_wrap(~season) +
    scale_x_discrete(labels=seq(2,nrow(sm),1)) +
    labs(x='sample size',y='standard deviation of lapse rate across years',title=paste0(tvar,' lapse rates for samples based on ', paste(vars,collapse=', ')))
  ggsave(paste0(figdir,tvar,'_similar_station_sd_vs_samplesize.jpeg'))
  # standard deviation (across years) of lapse rates increases with decreasing sample size
  
  plot(2:nrow(sm),elevrng)
  
  out <- lapse_tab %>% group_by(season,num) %>% summarise(mlapse=mean(lapse),sdl = sd(lapse))
  out %>% group_by(season) %>% summarise(sdcor = cor(sdl,elevrng)^2) # cor of the elevation range with the lapse rate standard deviation
  out %>% group_by(season) %>% summarise(mcor = cor(mlapse,elevrng)^2) # cor of the elvation range with the mean lapse rate
  
  ggplot(lapse_tab) +geom_line(aes(x=year,y=lapse_rmse,group=num,col=num)) + 
    labs(y='RMSE of lapse rate regression', title=paste0(tvar,' lapse rates RMSEs for samples based on ', paste(vars,collapse=', '))) +
    scale_color_manual(name='sample size',labels=seq(2,nrow(sm),1),values=c(rainbow(length(samp_size)-1),'black')) +
    facet_wrap(~season)
  ggsave(paste0(figdir,tvar,'_similar_station_RMSE_vs_samplesize.jpeg'))
  
  lapse_tab %>% group_by(season,num) %>% summarise(mrmse = mean(lapse_rmse)) %>%
    ggplot() +geom_point(aes(x=num,y=mrmse,col=season)) + 
    labs(y='mean RMSE of lapse rate regression', x='sample size', title=paste0(tvar,' lapse rates RMSEs for samples based on ', paste(vars,collapse=', ')))
  ggsave(paste0(figdir,tvar,'_similar_station_mRMSE_vs_samplesize.jpeg'))
  
  lapse_tab %>% group_by(season,num) %>% 
    summarise(ml=mean(lapse), me=mean(as.numeric(levels(elevrng))[elevrng])) %>%
    ggplot() + geom_point(aes(x=me,y=ml,col=num,group=num)) + facet_wrap(~season) +
    labs(x='Elevation Range (m)', y='Mean Lapse Rate', title=paste0(tvar,' lapse rates RMSEs for samples based on ', paste(vars,collapse=', ')))
  ggsave(paste0(figdir,tvar,'_similar_station_elevrng_vs_lapse.jpeg'))
  
  return(lapse_tab)
  
}


