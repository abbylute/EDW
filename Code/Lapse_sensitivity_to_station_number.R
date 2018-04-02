# How does lapse rate value vary depending on number of stations used to calculate it?

# set up some variables:
season_list <- unique(seasons$season)
years <- unique(season_means$year)
stations <- data.frame(unique(season_means$station_id))
snum <- c(3,5,10,15,20,25,30) # station sample sizes to evaluate
its <- 1000 # number of iterations

# set up the output table:
lout <- data.frame(matrix(ncol=3,nrow=(4*length(snum)*its)))
names(lout) <- c('season','station_num','lapse')
lout$season <- rep(season_list, each=(nrow(lout)/4))
lout$station_num <- rep(rep(snum, each=its), 4)

for (ss in 1:4){ # for each season
  sm <- season_means %>% filter(season==season_list[ss])
  #for (yy in 1:length(years)){

  for (nn in 1:length(snum)){
    for (ii in 1:its){
      xx <- which(lout$season==season_list[ss] & lout$station_num==snum[nn])[ii]
      
      tt <- sample_n(stations, size=snum[nn])
      tab <- sm %>% filter(station_id %in% tt[,1])
      lout$lapse[xx] <- summary(lm(data=tab, formula=smtemp ~ elev ))$coefficients[2]*1000
    } # end iterations
    
  } # end snum
  
  #} # end years
} # end ss

lout$station_num <- factor(lout$station_num,unique(lout$station_num))

ggplot(lout) +
  geom_boxplot(aes(x=station_num,y=lapse)) + facet_wrap(~season) +
  labs(x='number of stations',y='lapse rate (C/km)', 
       title='Tmin lapse rates',
       subtitle= 'calculated from 1000 random samples of a 49 station dataset')
ggsave(paste0(figdir,'Minimum_stations.jpeg'))
