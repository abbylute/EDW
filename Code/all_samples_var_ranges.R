all_samples_var_ranges <- function(sampcol, meta., meta_sc., vars.){
  library(tidyverse)
  
  samp_vars <- (matrix(nrow=1, ncol=(length(vars.)+6))) # set up an table with the ranges of each variable for each sample
  samp_var_names <- c('num_stations',vars.,'Annual','Spring','Summer','Fall','Winter')
  samp_vars[1] <- length(sampcol)
  cc <- which(names(meta.) %in% vars.)
  
  tt <- meta_sc.[which(meta_sc.$station_id %in% sampcol),]
  for (vv in 1:length(vars.)){
    samp_vars[(vv+1)] <- diff(range(tt[,cc[vv]])) # because I'm using the range function, there's no need to use the abs value, it does this (indirectly) automatically
  }
  df_small <- meta. %>% filter(station_id %in% sampcol) %>% group_by(season, station_id, elev) %>% summarise(mt = mean(tmax))
  ls <- df_small %>% group_by(season) %>% summarise(lapse = summary(lm(mt~elev))$coefficients[2]*1000)
  samp_vars[((ncol(samp_vars)-4):ncol(samp_vars))] <- ls$lapse[match(samp_var_names[((ncol(samp_vars)-4):ncol(samp_vars))], as.character(ls$season))]
  return(samp_vars)
}