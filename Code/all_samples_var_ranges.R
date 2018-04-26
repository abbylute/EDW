all_samples_var_ranges2 <- function(sampcol, meta., meta_sc., vars.){
  library(tidyverse)
  
  samp_vars <- (matrix(nrow=1, ncol=(length(vars.)+1+5+5))) # set up an table with the ranges of each variable for each sample
  samp_var_names <- c('num_stations',vars.,'Annual','Spring','Summer','Fall','Winter','Annual_rmse','Spring_rmse','Summer_rmse','Fall_rmse','Winter_rmse')
  samp_vars[1] <- length(sampcol) # this is the number of stations
  cc <- which(names(meta.) %in% vars.)
  
  tt <- meta_sc.[which(meta_sc.$station_id %in% sampcol),]
  for (vv in 1:length(vars.)){
    samp_vars[(vv+1)] <- diff(range(tt[,cc[vv]])) # because I'm using the range function, there's no need to use the abs value, it does this (indirectly) automatically
  }
  df_small <- meta. %>% filter(station_id %in% sampcol) %>% group_by(season, station_id, elev) %>% summarise(mt = mean(tmax))
  ls <- df_small %>% group_by(season) %>% summarise(lapse = summary(lm(mt~elev))$coefficients[2]*1000,
                                                    rmse = sqrt(sum(summary(lm(mt~elev))$residuals^2)/length(sampcol)))
  samp_vars[(length(vars)+2):(length(vars)+6)] <- ls$lapse[match(samp_var_names[(length(vars)+2):(length(vars)+6)], as.character(ls$season))]
  
  samp_vars[(length(vars)+7):(length(vars)+11)] <- ls$rmse[match(samp_var_names[(length(vars)+2):(length(vars)+6)], as.character(ls$season))]
  
  
  
  return(samp_vars)
}

