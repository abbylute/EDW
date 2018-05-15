all_samples_var_ranges <- function(sampcol, meta., meta_sc., vars.){
  library(tidyverse)
  
  samp_vars <- (matrix(nrow=(length(unique(meta.$season))*length(vars.)), ncol=6)) # set up an table with the ranges of each variable for each sample
  samp_var_names <- c('num_stations','season','lapse','rmse','var','var_range')
  samp_vars[,1] <- length(sampcol) # this is the number of stations
  samp_vars[,2] <- rep(unique(as.character(meta.$season)), each=length(vars.))
  samp_vars[,5] <- rep(vars., length(unique(meta.$season)))
  cc <- which(names(meta.) %in% vars.)
  
  tt <- meta_sc.[which(meta_sc.$station_id %in% sampcol),]
  var_rngs <- numeric(length(vars.))
  for (vv in 1:length(vars.)){
    var_rngs[vv] <- diff(range(tt[,cc[vv]])) # because I'm using the range function, there's no need to use the abs value, it does this (indirectly) automatically
  }
  samp_vars[,6] <-rep(var_rngs, length(unique(meta.$season)))
  
  df_small <- meta. %>% filter(station_id %in% sampcol) %>% group_by(season, station_id, elev) %>% summarise(mt = mean(temp))
  ls <- df_small %>% group_by(season) %>% summarise(lapse = summary(lm(mt~elev))$coefficients[2]*1000,
                                                    rmse = sqrt(sum(summary(lm(mt~elev))$residuals^2)/length(sampcol)))
  samp_vars[,3] <- rep(ls$lapse[match(unique(as.character(meta.$season)), as.character(ls$season))], each=length(vars.))
  samp_vars[,4] <- rep(ls$rmse[match(unique(as.character(meta.$season)), as.character(ls$season))], each=length(vars.))
  
  
  
  return(samp_vars)
}

