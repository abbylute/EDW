# Multiobjective Optimization
# get df and cortab from Local_seasonal_lapse_analysis.R
library(parallel)
library(tictoc)
source('Code/all_samples_var_ranges.R')

tvar <- 'tmax'
meta <- df[which(!is.na(df$tmax)),]
total_stations <- length(unique(meta$station_id))
vars <- c('elev','srad','m200','m500','m1000','m5000')
opt <- c('max',rep('min',5))
seas <- unique(df$season)
# scale: if scale=T (default), then normalize the vars
scale=T

# how many possible combinations are there for each sample size?
choose(total_stations,c(2:total_stations))
sum(choose(total_stations, c(2:total_stations)))

# get scaled variables:
cc <- which(names(meta) %in% vars)
meta_sc <- meta
if (scale==T){
  meta_sc[,cc] <- scale(meta[,cc])
}

# set up of the dataframe to contain all of the sample ranges
samp_df <- (matrix(nrow=0,ncol=(length(vars)+6)))

# Loop through sample sizes:
#for (num_stations in 2:total_stations){
num_stations <- 3
  
  samp <- combn(unique(meta_sc$station_id),num_stations) # compute all possible combinations of this number of stations
  samp <- cbind(as.character(samp[1,]),as.character(samp[2,]))
  # may have to trim sample sets because they get very large for more than ~20 stations:
    #samp <- samp[,sample(c(1:ncol(samp)),pmin(ncol(samp),10000), replace=F),drop=F]
  
  
  cl <- makeCluster(7)
  clusterExport(cl, c("all_samples_var_ranges","meta","meta_sc","vars"))
  samp_vars <- parApply(cl,samp, 1, function(i) all_samples_var_ranges(sampcol=i,meta,meta_sc,vars))
  stopCluster(cl)
  samp_vars <- t(samp_vars)

  samp_df <- rbind(samp_df,samp_vars) # combine the tables for each num_stations value
    
#} # end num_stations
names(samp_df) <- c('num_stations',vars,'Annual','Spring','Summer','Fall','Winter')
    
    
samp_dflong <- samp_df %>% gather('season', 'lapse', Annual:Winter)
samp_dflong <- samp_dflong %>% gather('var','var_range', elev:m5000)
cors <- cortab %>% filter(tvar == 'tmax')
samp_dflong <- left_join(samp_dflong,cors, by=c('season','var'))
#samp_dflong$var_range[which(samp_dflong$var=='elev')] <- 

# rescale elevation ranges so that larger ranges get smaller numbers (like all the other variables)
evr <- samp_dflong$var_range[which(samp_dflong$var == 'elev')] # elevation ranges
ee <- ecdf(evr)
plot(ee) # the original distribution

  # a rescale option:
  devr <- max(evr)
  ee <- ecdf(devr-evr)
  plot(ee) # the rescaled distribution

samp_dflong$var_range[which(samp_dflong$var == 'elev')] <- devr-evr

# weight the ranges by the variable's correlation with seasonal temperature
# old way, caused unimportant variables to get large values, making it potentially harder to find the best samples
  samp_wt <- samp_dflong %>% mutate(var_range_wt = var_range*abs(1/cor)) # causes ranges to be larger when correlations are lower
  samp_wt_wide <- samp_wt %>% select(-tvar,-cor,-var_range) %>% spread('var','var_range_wt')
  samp_wt_wide <- samp_wt_wide %>% mutate(vsum = sqrt(elev^2 + m1000^2 + m200^2 + m500^2 + m5000^2 + srad^2))
  ggplot(samp_wt_wide) + geom_point(aes(x=vsum, y=lapse)) + facet_wrap(~season) + 
    labs(x='sum of weighted variable ranges',y='lapse rate', title='Tmax inverse weighted multiobjective optimization (n=5)',
         subtitle = 'variables include elevation, srad, and all TPIs')
  ggsave(paste0(figdir,'tmax_multi_opt_inversewt_n5.jpeg'))
  samp_wt_wide %>% group_by(season) %>% filter(vsum == min(vsum))
  

  # weight the ranges by the variable's correlation with seasonal temperature
  # second way, caused unimportant variables to get smaller values, making it maybe easier to ignore these in trying to find best sample
  samp_wt1 <- samp_dflong %>% mutate(var_range_wt = var_range*abs(cor)) # causes ranges to be larger when correlations are lower
  samp_wt_wide1 <- samp_wt1 %>% select(-tvar,-cor,-var_range) %>% spread('var','var_range_wt')
  samp_wt_wide1 <- samp_wt_wide1 %>% mutate(vsum = sqrt(elev^2 + m1000^2 + m200^2 + m500^2 + m5000^2 + srad^2))
  ggplot(samp_wt_wide1) + geom_point(aes(x=vsum, y=lapse)) + facet_wrap(~season) + 
    labs(x='sum of weighted variable ranges',y='lapse rate', title='Tmax weighted multiobjective optimization (n=5)',
         subtitle = 'variables include elevation, srad, and all TPIs')
  ggsave(paste0(figdir,'tmax_multi_opt_wt_n5.jpeg'))
  samp_wt_wide1 %>% group_by(season) %>% filter(vsum == min(vsum)) %>% summarise(mlapse=mean(lapse))
  
    
  # weight the ranges by the variable's correlation with seasonal temperature
  # new way, 
  samp_wt_wide2 <- samp_dflong %>% select(-tvar,-cor) %>% spread('var','var_range')
  outtab <- data.frame(matrix(ncol=(ncol(samp_wt_wide2)+1),nrow=0)); names(outtab) <- c(names(samp_wt_wide2),'vsum')
  for (ss in 1:length(seas)){
    tab <- samp_wt_wide2 %>% filter(season == seas[ss])
    corstrim <- cors %>% filter(season == seas[ss])
    tab <- tab %>% mutate(vsum = sqrt( (elev^2)*abs(corstrim$cor[which(corstrim$var=='elev')]) +
                                (srad^2)*abs(corstrim$cor[which(corstrim$var=='srad')]) +
                                (m200^2)*abs(corstrim$cor[which(corstrim$var=='m200')]) +
                                (m500^2)*abs(corstrim$cor[which(corstrim$var=='m500')]) +
                                (m1000^2)*abs(corstrim$cor[which(corstrim$var=='m1000')]) +
                                (m5000^2)*abs(corstrim$cor[which(corstrim$var=='m5000')])))
    outtab <- rbind(outtab,tab)
  }
  ggplot(outtab) + geom_point(aes(x=vsum, y=lapse)) + facet_wrap(~season) + 
    labs(x='sum of weighted variable ranges',y='lapse rate', title='Tmax later weighted multiobjective optimization (n=5)',
         subtitle = 'variables include elevation, srad, and all TPIs')
  ggsave(paste0(figdir,'tmax_multi_opt_latewt_n5.jpeg'))
  outtab %>% group_by(season) %>% filter(vsum == min(vsum)) %>% summarise(mlapse=mean(lapse))
  
  
  
  

# same but without the weighting:
samp_unwt <- samp_dflong %>% select(-tvar,-cor) %>% spread('var','var_range')
samp_unwt <- samp_unwt %>% mutate(vsum = sqrt(elev^2 + m1000^2 + m200^2 + m500^2 + m5000^2 + srad^2))
ggplot(samp_unwt) + geom_point(aes(x=vsum, y=lapse)) + facet_wrap(~season) + 
  labs(x='sum of unweighted variable ranges',y='lapse rate', title='Tmax unweighted multiobjective optimization (n=5)',
       subtitle = 'variables include elevation, srad, and all TPIs')
ggsave(paste0(figdir,'tmax_multi_opt_unwt_n5.jpeg'))
samp_unwt %>% group_by(season) %>% filter(vsum == min(vsum))

# also look at how each variable's ranges compares with the calculated lapse rates:
ggplot(samp_unwt) + geom_point(aes(x=elev, y=lapse)) + facet_wrap(~season) # strong funnel
ggsave(paste0(figdir,'tmax_multi_opt_unwt_elev_n5.jpeg'))
ggplot(samp_unwt) + geom_point(aes(x=m1000, y=lapse)) + facet_wrap(~season) # no clear pattern
ggsave(paste0(figdir,'tmax_multi_opt_unwt_m1000_n5.jpeg'))
ggplot(samp_unwt) + geom_point(aes(x=m200, y=lapse)) + facet_wrap(~season) # no clear pattern
ggsave(paste0(figdir,'tmax_multi_opt_unwt_m200_n5.jpeg'))
ggplot(samp_unwt) + geom_point(aes(x=m500, y=lapse)) + facet_wrap(~season) # no clear pattern
ggsave(paste0(figdir,'tmax_multi_opt_unwt_m500_n5.jpeg'))
ggplot(samp_unwt) + geom_point(aes(x=m5000, y=lapse)) + facet_wrap(~season) # no clear pattern
ggsave(paste0(figdir,'tmax_multi_opt_unwt_m5000_n5.jpeg'))
ggplot(samp_unwt) + geom_point(aes(x=srad, y=lapse)) + facet_wrap(~season) # no clear pattern
ggsave(paste0(figdir,'tmax_multi_opt_unwt_srad_n5.jpeg'))
# elevation is the only variable that helps collapse the range of potential lapse rate values.




dat <- samp_wt_wide %>% filter(season=='Annual')
plot(dat$vsum, dat$lapse)

library(colorRamps)
ggplot(dat, aes(x=elev, y=m500, z=srad)) + geom_point(aes(col=lapse)) +
  scale_color_gradientn(name='lapse rate\n(C/km)',colors=matlab.like(50))


library(plotly)
colors <- matlab.like(100)
p <- plot_ly(dat, x = ~elev, y = ~m500, z = ~srad, color = ~lapse, colors = colors)
p







