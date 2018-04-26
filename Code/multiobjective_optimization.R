# Multiobjective Optimization
# get df and cortab from Local_seasonal_lapse_analysis.R first
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
# took 50 minutes to do a sample size of 20 (2:20)
tic('all samples')
for (num_stations in 2:total_stations){
#num_stations <- 25
  
  samp <- t(combn(as.character(unique(meta_sc$station_id)),num_stations)) # compute all possible combinations of this number of stations
   # may have to trim sample sets because they get very large for more than ~20 stations:
    #samp <- samp[,sample(c(1:ncol(samp)),pmin(ncol(samp),10000), replace=F),drop=F]
  
  
  cl <- makeCluster(7)
  clusterExport(cl, c("all_samples_var_ranges","meta","meta_sc","vars"))
  tic('go')
  print(num_stations)
  samp_vars <- parApply(cl,samp, 1, function(i) all_samples_var_ranges(sampcol=i,meta,meta_sc,vars))
  toc()
  stopCluster(cl)
  samp_vars <- t(samp_vars)

  samp_df <- rbind(samp_df,samp_vars) # combine the tables for each num_stations value
    
} # end num_stations
toc()
samp_df <- data.frame(samp_df)
names(samp_df) <- c('num_stations',vars,'Annual','Spring','Summer','Fall','Winter')
    write.table(samp_df, paste0(figdir,gsub(' ','',regname),'_allsamples_data.csv'))
    
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

### PLOTTING SUMMED WEIGHTS ###
wt_opt <-   2 # specify a wt_scheme
wt_schemes <- c('unwt','latewt','inversewt','wt')
wt_name <- c('unweighted','late weighted','inverse weighted','weighted')

nn <- length(unique(samp_dflong$num_stations))
  if (wt_opt ==1){
    samp_wt_wide <- samp_dflong %>% select(-tvar,-cor) %>% spread('var','var_range')
    samp_wt_wide <- samp_wt_wide %>% mutate(vsum = sqrt(elev^2 + m1000^2 + m200^2 + m500^2 + m5000^2 + srad^2))
    
  } else if (wt_opt==2){
    samp_wt <- samp_dflong %>% select(-tvar,-cor) %>% spread('var','var_range')
    samp_wt_wide <- data.frame(matrix(ncol=(ncol(samp_wt)+1),nrow=0)); names(samp_wt_wide) <- c(names(samp_wt),'vsum')
    for (ss in 1:length(seas)){
      tab <- samp_wt %>% filter(season == seas[ss])
      corstrim <- cors %>% filter(season == seas[ss])
      tab <- tab %>% mutate(vsum = sqrt((elev^2)*abs(corstrim$cor[which(corstrim$var=='elev')]) + (srad^2)*abs(corstrim$cor[which(corstrim$var=='srad')]) +
                                        (m200^2)*abs(corstrim$cor[which(corstrim$var=='m200')]) + (m500^2)*abs(corstrim$cor[which(corstrim$var=='m500')]) +
                                        (m1000^2)*abs(corstrim$cor[which(corstrim$var=='m1000')]) + (m5000^2)*abs(corstrim$cor[which(corstrim$var=='m5000')])))
      samp_wt_wide <- rbind(samp_wt_wide,tab)
    }
    
  } else if (wt_opt==3){
    samp_wt <- samp_dflong %>% mutate(var_range_wt = var_range*abs(1/cor))  # causes ranges to be larger when correlations are lower
    samp_wt_wide <- samp_wt  %>% select(-tvar,-cor,-var_range) %>% spread('var','var_range_wt')
    samp_wt_wide <- samp_wt_wide %>% mutate(vsum = sqrt(elev^2 + m1000^2 + m200^2 + m500^2 + m5000^2 + srad^2))
    
  } else if (wt_opt==4){
    samp_wt <- samp_dflong %>% mutate(var_range_wt = var_range*abs(cor)) # causes ranges to be larger when correlations are lower
    samp_wt_wide <- samp_wt %>% select(-tvar,-cor,-var_range) %>% spread('var','var_range_wt')
    samp_wt_wide <- samp_wt_wide %>% mutate(vsum = sqrt(elev^2 + m1000^2 + m200^2 + m500^2 + m5000^2 + srad^2))
    
  }
  
  # plot it:
  if(nn ==1){
    ggplot(samp_wt_wide) + geom_point(aes(x=vsum, y=lapse)) + facet_wrap(~season) + 
      labs(x='sum of weighted variable ranges',y='lapse rate', title=paste0('Tmax ',wt_name[wt_opt],' multiobjective optimization (n=',num_stations,')'),
           subtitle = 'variables include elevation, srad, and all TPIs')
    #ggsave(paste0(figdir,'tmax_multi_opt_',wt_schemes[wt_opt],'_n',num_stations,'.jpeg'))
    samp_wt_wide %>% group_by(season) %>% filter(vsum == min(vsum))
    
  } else { # for multiple sample sizes
    for (ss in 1:length(seas)){
      sname <- as.character(seas[ss])
      ggplot(samp_wt_wide[which(samp_wt_wide$season==sname),]) + geom_point(aes(x=vsum, y=lapse)) + facet_wrap(~num_stations) + 
        lims(y=c(-30,20)) +
        labs(x='sum of weighted variable ranges',y='lapse rate', title=paste0(sname,' Tmax ',wt_name[wt_opt],' multiobjective optimization'),
             subtitle = 'variables include elevation, srad, and all TPIs')
      ggsave(paste0(figdir,sname,'_tmax_multi_opt_',wt_schemes[wt_opt],'_allsamples.jpeg'))
      
    } #seasons
  } # plotting
####### end plotting summed weights ########




### PLOTTING INDIVIDUAL VARIABLE WEIGHTS ###
nn <- length(unique(samp_dflong$num_stations))
# simply use unweighted:
samp_wt_wide <- samp_dflong %>% select(-tvar,-cor) %>% spread('var','var_range')
samp_wt_wide <- samp_wt_wide %>% mutate(vsum = sqrt(elev^2 + m1000^2 + m200^2 + m500^2 + m5000^2 + srad^2))

if (nn==1){
  for (vv in 1:length(vars)){
    ggplot(samp_wt_wide) + geom_point(aes_string(x=vars[vv], y='lapse')) + 
      lims(y=c(-30,20)) + facet_wrap(~season) +
      labs(x=paste0(vars[vv],' range'),y='lapse rate', title=paste0('Tmax ',wt_name[wt_opt],' multiobjective optimization'))
    ggsave(paste0(figdir,'tmax_multi_opt_unwt_',vars[vv],'_n',num_stations,'.jpeg'))
   
  } # end variables
} else { # for multiple sample sizes
  for (vv in 1:length(vars)){
    for (ss in 1:length(seas)){
      sname <- as.character(seas[ss])
      sampy <- samp_wt_wide[which(samp_wt_wide$season==seas[ss]),]
      names(sampy)[which(names(sampy)==vars[vv])] <- 'vname'
      minlapse <- sampy %>% group_by(num_stations) %>% filter(vname == min(vname)) %>% select(num_stations,lapse)
      minlapse$lapse <- round(minlapse$lapse,2)
      minlapse$x <- rep(min(sampy$vname)+.2*(diff(range(sampy$vname))),nrow(minlapse))
      minlapse$y <- rep(18, nrow(minlapse))
      
      ggplot(sampy) + geom_point(aes(x=vname, y=lapse)) + geom_text(data=minlapse, aes(x=x,y=y,label=lapse)) +
        lims(y=c(-30,20)) + facet_wrap(~num_stations) +
        labs(x=paste0(vars[vv],' range'),y='lapse rate', title=paste0(sname,' tmax ',wt_name[wt_opt],' multiobjective optimization'))
      ggsave(paste0(figdir,sname,'_tmax_multi_opt_unwt_',vars[vv],'_allsamples.jpeg'))
    
    } # end seasons
  } # end variables
} # end plotting
rm(sampy);gc()
####### end plotting individual variable weights ########





dat <- samp_wt_wide %>% filter(season=='Annual')
plot(dat$vsum, dat$lapse)

library(colorRamps)
ggplot(dat, aes(x=elev, y=m500, z=srad)) + geom_point(aes(col=lapse)) +
  scale_color_gradientn(name='lapse rate\n(C/km)',colors=matlab.like(50))


library(plotly)
colors <- matlab.like(100)
p <- plot_ly(dat, x = ~elev, y = ~m500, z = ~srad, color = ~lapse, colors = colors)
p







