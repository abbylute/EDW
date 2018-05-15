# Multiobjective Optimization
# get df and cortab from Local_seasonal_lapse_analysis.R first
library(parallel)
library(tictoc)
library(colorRamps)
source('Code/all_samples_var_ranges.R')

tvar <- 'tmax'
names(df)[which(names(df)==tvar)] <- 'temp'
meta <- df %>% filter(!is.na(tvar))
total_stations <- length(unique(meta$station_id))
vars <- c('elev','srad','m200','m500','m1000','m5000','windex_1000')
opt <- c('max',rep('min',(length(vars)-1)))
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
samp_df <- (matrix(nrow=0,ncol=6))

# Loop through sample sizes:
# took 50 minutes to do a sample size of 20 (2:20)
tic('all samples')
for (num_stations in 2:total_stations){
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
  
  # reshape sampvars:
  rowcount <- length(unique(meta$season))*length(vars)
  outtab <- matrix(nrow=(rowcount*ncol(samp_vars)), ncol=ncol(samp_df))
  for (dd in 1:ncol(samp_df)){
    theserows <- ((dd-1)*rowcount+1):(dd*rowcount)
    outtab[,dd] <- matrix(samp_vars[theserows,],ncol=1,byrow=T)
    
  }
  
  samp_df <- rbind(samp_df,outtab) # combine the tables for each num_stations value
  
} # end num_stations
toc()
samp_df <- data.frame(samp_df, stringsAsFactors=F)

names(samp_df) <- c('num_stations','season','lapse','rmse','var','var_range')
class(samp_df[,1]) <- 'numeric'
class(samp_df[,3]) <- 'numeric'
class(samp_df[,4]) <- 'numeric'
class(samp_df[,6]) <- 'numeric'

#write.table(samp_df, paste0(figdir,gsub(' ','',regname),'_allsamples_data_wrmse.csv'))
#samp_df <- read.table(paste0(figdir,gsub(' ','',regname),'_allsamples_data_wrmse.csv'))
class(samp_df[,1]) <- 'numeric'
cors <- cortab %>% filter(tvar == 'tmax')
samp_df <- left_join(samp_df,cors, by=c('season','var')) %>% dplyr::select(-tvar)

# rescale elevation ranges so that larger ranges get smaller numbers (like all the other variables)
evr <- samp_df$var_range[which(samp_df$var == 'elev')] # elevation ranges
ee <- ecdf(evr)
plot(ee) # the original distribution

# a rescale option:
devr <- max(evr)
ee <- ecdf(devr-evr)
plot(ee) # the rescaled distribution

samp_df$var_range[which(samp_df$var == 'elev')] <- devr-evr

### PLOTTING SUMMED WEIGHTS ###
wt_opt <-   1 # specify a wt_scheme
wt_schemes <- c('unwt','latewt','inversewt','wt')
wt_name <- c('unweighted','late weighted','inverse weighted','weighted')

nn <- length(unique(samp_df$num_stations))
if (wt_opt ==1){
  samp_wt_wide <- samp_df %>% dplyr::select(-cor) %>% spread('var','var_range')
  samp_wt_wide <- samp_wt_wide %>% mutate(vsum = sqrt(elev^2 + m1000^2 + m200^2 + m500^2 + m5000^2 + srad^2 + windex_1000^2))
  
} else if (wt_opt==2){
  samp_wt <- samp_df %>% dplyr::select(-cor) %>% spread('var','var_range')
  samp_wt_wide <- data.frame(matrix(ncol=(ncol(samp_wt)+1),nrow=0)); names(samp_wt_wide) <- c(names(samp_wt),'vsum')
  for (ss in 1:length(seas)){
    tab <- samp_wt %>% filter(season == seas[ss])
    corstrim <- cors %>% filter(season == seas[ss])
    tab <- tab %>% mutate(vsum = sqrt((elev^2)*abs(corstrim$cor[which(corstrim$var=='elev')]) + (srad^2)*abs(corstrim$cor[which(corstrim$var=='srad')]) +
                                        (m200^2)*abs(corstrim$cor[which(corstrim$var=='m200')]) + (m500^2)*abs(corstrim$cor[which(corstrim$var=='m500')]) +
                                        (m1000^2)*abs(corstrim$cor[which(corstrim$var=='m1000')]) + (m5000^2)*abs(corstrim$cor[which(corstrim$var=='m5000')]) +
                                        (windex_1000^2)*abs(corstrim$cor[which(corstrim$var=='windex_1000')])))
    samp_wt_wide <- rbind(samp_wt_wide,tab)
  }
  
} else if (wt_opt==3){
  samp_wt <- samp_df %>% mutate(var_range_wt = var_range*abs(1/cor))  # causes ranges to be larger when correlations are lower
  samp_wt_wide <- samp_wt  %>% dplyr::select(-cor,-var_range) %>% spread('var','var_range_wt')
  samp_wt_wide <- samp_wt_wide %>% mutate(vsum = sqrt(elev^2 + m1000^2 + m200^2 + m500^2 + m5000^2 + srad^2 + windex_1000^2))
  
} else if (wt_opt==4){
  samp_wt <- samp_df %>% mutate(var_range_wt = var_range*abs(cor)) # causes ranges to be larger when correlations are lower
  samp_wt_wide <- samp_wt %>% dplyr::select(-cor,-var_range) %>% spread('var','var_range_wt')
  samp_wt_wide <- samp_wt_wide %>% mutate(vsum = sqrt(elev^2 + m1000^2 + m200^2 + m500^2 + m5000^2 + srad^2 + windex_1000^2))
  
}

# plot it:
if(nn ==1){
  ggplot(samp_wt_wide) + geom_point(aes(x=vsum, y=lapse, col=rmse)) + facet_wrap(~season) + 
    labs(x='sum of weighted variable ranges',y='lapse rate', title=paste0('Tmax ',wt_name[wt_opt],' multiobjective optimization (n=',num_stations,')'),
         subtitle = 'variables include elevation, srad, and all TPIs') + scale_color_gradientn(colours=matlab.like(50))
  #ggsave(paste0(figdir,'tmax_multi_opt_',wt_schemes[wt_opt],'_n',num_stations,'.jpeg'))
  samp_wt_wide %>% group_by(season) %>% filter(vsum == min(vsum))
  
} else { # for multiple sample sizes
  for (ss in 1:length(seas)){
    sname <- as.character(seas[ss])
    print(sname)
    jpeg(filename=paste0(figdir,sname,'_tmax_multi_opt_',wt_schemes[wt_opt],'_allsamples.jpeg'),width=10,height=10,unit='in', res=500)
    gg <- ggplot(samp_wt_wide[which(samp_wt_wide$season==sname),]) + geom_point(aes(x=vsum, y=lapse, col=rmse),size=.7) + facet_wrap(~num_stations) + 
      lims(y=c(-30,20)) + scale_color_gradientn(colours=matlab.like(50)) + 
      theme(panel.background = element_rect(fill = 'white'), panel.grid.major = element_line(colour = "lightgrey"),
            panel.grid.minor=element_line(colour="lightgrey")) + #, panel.border = element_rect(colour = "lightgrey", fill=NA, size=1)) +
      labs(x='sum of weighted variable ranges',y='lapse rate', title=paste0(sname,' Tmax ',wt_name[wt_opt],' multiobjective optimization'),
           subtitle = 'variables include elevation, srad, and all TPIs')
    print(gg)
    dev.off()
    #ggsave(paste0(figdir,sname,'_tmax_multi_opt_',wt_schemes[wt_opt],'_allsamples.jpeg'))
   
    #plot rmses
     gg <- ggplot(samp_wt_wide[which(samp_wt_wide$season==sname),]) + geom_point(aes(x=vsum, y=rmse, col=lapse),size=.7) + 
      facet_wrap(~num_stations) + 
      #lims(y=c(-30,20)) + 
      scale_color_gradientn(colours=matlab.like(50),limits=c(-15,5)) + 
      theme(panel.background = element_rect(fill = 'white'), panel.grid.major = element_line(colour = "lightgrey"),
            panel.grid.minor=element_line(colour="lightgrey")) + #, panel.border = element_rect(colour = "lightgrey", fill=NA, size=1)) +
      labs(x='sum of weighted variable ranges',y='rmse', title=paste0(sname,' Tmax ',wt_name[wt_opt],' multiobjective optimization'),
           subtitle = 'variables include elevation, srad, and all TPIs')
      print(gg)
      #ggsave(paste0(figdir,sname,'_tmax_multi_opt_',wt_schemes[wt_opt],'_allsamples_rmse.jpeg'))
      
    
  } #seasons
} # plotting
####### end plotting summed weights ########


#### Experiment with plotting pdfs of lapse rates for categories of vsum values: ####
  sea <- "Annual"  
  samp_cat <- samp_wt_wide %>% filter(season==sea) %>% mutate(vsumcat = ceiling(vsum), rmsecat = ceiling(rmse)) %>% 
    dplyr::select(num_stations,lapse,vsumcat,rmsecat) %>% ungroup()
  samp_cat$vsumcat <- factor(samp_cat$vsumcat, rev(sort(unique(samp_cat$vsumcat))))
  samp_cat$rmsecat <- factor(samp_cat$rmsecat, rev(sort(unique(samp_cat$rmsecat))))
  samp_cat$lapse <- pmin(samp_cat$lapse,10); samp_cat$lapse <- pmax(samp_cat$lapse,-20)
  
  # vsum plots:
  # version with free scales:
  ggplot(samp_cat) +
    geom_density(aes(x=lapse, fill=vsumcat), alpha=.6) + facet_wrap(~num_stations, scales="free") +
    labs(x='lapse rate (C/km)', title=paste0(sea,' ',wt_name[wt_opt],' pdfs based on sum of ranges')) +
    scale_fill_discrete(name='sum of ranges\ncategory')
  ggsave(paste0(figdir,sea,'_',wt_schemes[wt_opt],'_lapse_density_by_vsum_freescale.jpeg'))
  # version with fixed scales:
  ggplot(samp_cat) +
    geom_density(aes(x=lapse, fill=vsumcat), alpha=.6) + facet_wrap(~num_stations) +
    labs(x='lapse rate (C/km)', title=paste0(sea,' ',wt_name[wt_opt],' pdfs based on sum of ranges')) +
    scale_fill_discrete(name='sum of ranges\ncategory')
  ggsave(paste0(figdir,sea,'_',wt_schemes[wt_opt],'_lapse_density_by_vsum_fixedscale.jpeg'))

  # rmse plots:
  # version with free scales:
  ggplot(samp_cat) +
    geom_density(aes(x=lapse, fill=rmsecat), alpha=.6) + facet_wrap(~num_stations, scales="free") +
    labs(x='lapse rate (C/km)', title=paste0(sea,' ',wt_name[wt_opt],' pdfs based on rmse')) +
    scale_fill_discrete(name='rmse\ncategory')
  ggsave(paste0(figdir,sea,'_',wt_schemes[wt_opt],'_lapse_density_by_rmse_freescale.jpeg'))
  # version with fixed scales:
  ggplot(samp_cat) +
    geom_density(aes(x=lapse, fill=vsumcat), alpha=.6) + facet_wrap(~num_stations) +
    labs(x='lapse rate (C/km)', title=paste0(sea,' ',wt_name[wt_opt],' pdfs based on rmse')) +
    scale_fill_discrete(name='rmse\ncategory')
  ggsave(paste0(figdir,sea,'_',wt_schemes[wt_opt],'_lapse_density_by_rmse_fixedscale.jpeg'))
  
  #look at how rmse and vsum compare:
  ggplot() +
    geom_density(data=samp_cat[which(samp_cat$rmsecat==1),], aes(x=lapse, fill=rmsecat), fill='blue', alpha=.6) + 
    geom_density(data=samp_cat[which(samp_cat$vsumcat==1),], aes(x=lapse, fill=vsumcat), fill='red', alpha=.6) + 
    facet_wrap(~num_stations, scales="free") +
    labs(x='lapse rate (C/km)', title=paste0(sea,' ',wt_name[wt_opt],' pdfs based on rmse and sum of ranges'),
         subtitle='blue=lowest rmse category, red=lowest sum of ranges category') 
  ggsave(paste0(figdir,sea,'_',wt_schemes[wt_opt],'_lapse_density_by_rmse_and_vsum.jpeg'))
  
  

#### Plot lapse rates for lowest vsum, lowest rmse for each sample size and season. With % range ####
  prng <- .05
minlapse1 <- samp_wt_wide %>% group_by(num_stations,season) %>% filter(vsum == min(vsum)) %>% dplyr::select(num_stations,lapse,vsum,rmse) %>% 
  group_by(num_stations,season,vsum) %>% summarise(lapse=round(mean(lapse),2), rmse=round(mean(rmse),2)) # add summarise call in case there are multiple instances of the minimum value
minlapse2 <- samp_wt_wide %>% group_by(num_stations,season) %>% filter(rmse == min(rmse)) %>% dplyr::select(num_stations,lapse,vsum,rmse) %>%
  group_by(num_stations,season) %>% summarise(lapse=round(mean(lapse),2), rmse=round(mean(rmse),2)) %>% filter(!num_stations==2)
minlapse1r <- samp_wt_wide %>% group_by(num_stations,season) %>% mutate(ncount =round(n()*prng)) %>% 
  group_by(num_stations,season) %>% filter(vsum %in% sort(vsum)[1:ncount]) %>% group_by(num_stations,season) %>% 
  summarise(ymin=min(lapse), ymax=max(lapse))
minlapse2r <- samp_wt_wide %>% group_by(num_stations,season) %>% mutate(ncount=round(n()*prng)) %>%
  group_by(num_stations,season) %>% filter(rmse %in% sort(rmse)[1:ncount]) %>% 
  summarise(ymin=min(lapse), ymax=max(lapse)) %>% filter(!num_stations == 2)

ggplot(minlapse1) + geom_line(aes(x=num_stations,y=lapse, group=season)) + 
  geom_line(inherit.aes=F, data=minlapse2, aes(x=num_stations, y=lapse, group=season),col='red') +
  geom_ribbon(inherit.aes=F, data=minlapse1r, aes(x=num_stations, ymin=ymin, ymax=ymax, group=season), alpha=.5) +
  geom_ribbon(inherit.aes=F, data=minlapse2r, aes(x=num_stations, ymin=ymin, ymax=ymax, group=season), alpha=.5, fill='red') +
  facet_wrap(~season) +
  labs(title='Lapse rate of most similar station sample (Eastern Oregon)', 
       subtitle=(paste0('black line is lowest sum of ranges, red line is lowest rmse\nshaded areas represent the lowest ',prng*100,'% of sums/rmses'))) +
  coord_cartesian(ylim=c(-8.5,1))
ggsave(paste0(figdir,wt_schemes[wt_opt],'_season_lapse_for_mostsimilarsample.jpeg'))
rm(minlapse1,minlapse2,minlapse1r,minlapse2r); gc()

#ggplot(minlapse) + geom_point(aes(x=num_stations, y=vsum)) + labs(title=wt_name[wt_opt],y='station similarity',x='sample size') +
##  geom_point(aes(x=num_stations,y=rmse,col=season)) +lims(y=c(0,max(minlapse$rmse,minlapse$vsum))) +
#  scale_color_discrete(name='RMSE')
#ggsave(paste0(figdir,wt_schemes[wt_opt],'_similarity_v_samplesize.jpeg'))

minl <- samp_wt_wide %>% group_by(num_stations,season) %>% mutate(ncount =round(n()*prng)) %>% 
  group_by(num_stations,season) %>% filter(vsum %in% sort(vsum)[1:ncount]) %>% group_by(num_stations,season) %>% 
  summarise(ymin=min(lapse), ymax=max(lapse), yrng =ymax-ymin, mvsum=mean(vsum))

ggplot(minl) + geom_point(aes(x=mvsum, y=yrng, col=num_stations)) + facet_wrap(~season) +
  scale_color_gradientn(colors=matlab.like(19))


#### PLOTTING INDIVIDUAL VARIABLE WEIGHTS ###
nn <- length(unique(samp_df$num_stations))
# simply use unweighted:
samp_wt_wide <- samp_df %>% select(-cor) %>% spread('var','var_range')
samp_wt_wide <- samp_wt_wide %>% mutate(vsum = sqrt(elev^2 + m1000^2 + m200^2 + m500^2 + m5000^2 + srad^2))

if (nn==1){
  for (vv in 1:length(vars)){
    ggplot(samp_wt_wide) + geom_point(aes_string(x=vars[vv], y='lapse', col='rmse')) + 
      lims(y=c(-30,20)) + facet_wrap(~season) + scale_color_gradientn(colours=matlab.like(50)) +
      labs(x=paste0(vars[vv],' range'),y='lapse rate', title=paste0('Tmax ',wt_name[wt_opt],' multiobjective optimization'))
    ggsave(paste0(figdir,'tmax_multi_opt_unwt_',vars[vv],'_n',num_stations,'.jpeg'))
    
  } # end variables
} else { # for multiple sample sizes
  for (vv in 1:length(vars)){
    for (ss in 1:length(seas)){
      sname <- as.character(seas[ss])
      sampy <- samp_wt_wide[which(samp_wt_wide$season==seas[ss]),]
      names(sampy)[which(names(sampy)==vars[vv])] <- 'vname'
      minlapse <- sampy %>% group_by(num_stations) %>% filter(vname == min(vname)) %>% select(num_stations,lapse) %>% 
        group_by(num_stations) %>% summarise(lapse=mean(lapse)) # add summarise call in case there are multiple instances of the minimum value
      minlapse$lapse <- round(minlapse$lapse,2)
      minlapse$x <- rep(min(sampy$vname)+.2*(diff(range(sampy$vname))),nrow(minlapse))
      minlapse$y <- rep(18, nrow(minlapse))
      
      jpeg(filename=paste0(figdir,sname,'_tmax_multi_opt_unwt_',vars[vv],'_allsamples.jpeg'), width=10,height=10,unit='in', res=500)
      gg <-ggplot(sampy) + geom_point(aes(x=vname, y=lapse, col=rmse), size=.5) + geom_text(data=minlapse, aes(x=x,y=y,label=lapse)) +
        lims(y=c(-30,20)) + facet_wrap(~num_stations) + scale_color_gradientn(colours=matlab.like(50)) +
        theme(panel.background = element_rect(fill = 'white'), panel.grid.major = element_line(colour = "lightgrey"),
              panel.grid.minor=element_line(colour="lightgrey")) + #, panel.border = element_rect(colour = "lightgrey", fill=NA, size=1)) +
        labs(x=paste0(vars[vv],' range'),y='lapse rate', title=paste0(sname,' tmax ',wt_name[wt_opt],' multiobjective optimization, based on ',vars[vv]))
      print(gg)
      dev.off()
      #ggsave(paste0(figdir,sname,'_tmax_multi_opt_unwt_',vars[vv],'_allsamples.jpeg'))
      
    } # end seasons
  } # end variables
} # end plotting
rm(sampy);gc()
####### end plotting individual variable weights ########



#### Which samples are being selected? ####
# simply use unweighted:
#samp_wt_wide <- samp_df %>% select(-cor) %>% spread('var','var_range')
#samp_wt_wide <- samp_wt_wide %>% mutate(vsum = sqrt(elev^2 + m1000^2 + m200^2 + m500^2 + m5000^2 + srad^2))
mins <- samp_wt_wide %>% filter(season=='Fall') %>% group_by(num_stations) %>% summarise(vmin =which.min(vsum))

min_st <- data.frame(matrix(nrow=0,ncol=2))
for (nn in 2:total_stations){
  nns <- rep(nn,nn)
  station_id <- t(combn(as.character(unique(meta_sc$station_id)),nn))[mins$vmin[nn-1],]
  
  min_st <- rbind(min_st, cbind(nns,station_id))
}

stations <- data.frame(cbind(as.character(unique(meta$station_id)),1:20))
names(stations) <- c('station_id','id')

meta_trim <- meta %>% filter(season=='Annual') %>% group_by(station_id,lat,lon,elev,m200,m500,m1000,m5000,srad) %>% summarise(tmax=mean(tmax),tmin=mean(tmin))

min_st <- left_join(min_st,stations, by='station_id')
min_st <- left_join(min_st, meta_trim, by='station_id')
ggplot(min_st) + geom_point(aes(x=nns,y=id,col=elev)) + labs(x='sample size',y='station id') +
  scale_color_gradientn(colors=matlab.like(50))
# this is dropping the lowest valleys (m200) from most samples

ggplot() + annotate(geom='point', x=meta_trim$tmax,y=meta_trim$elev) +
geom_point(data=min_st, aes(x=tmax,y=elev),col='red') +facet_wrap(~nns) +
  labs(x='Annual Tmax (C)',y='elevation (m)',title='Eastern Oregon stations selected by multiobjective optimization',
       subtitle='red is a selected station for each sample size')
#ggsave(paste0(figdir,'stations_selected_by_multiobj_opt.jpeg'))

ggplot() + annotate(geom='point', x=meta_trim$tmax,y=meta_trim$m500) +
  geom_point(data=min_st, aes(x=tmax,y=m500),col='red') +facet_wrap(~nns) +
  labs(x='Annual Tmax (C)',y='TPI 500m',title='Eastern Oregon stations selected by multiobjective optimization',
       subtitle='red is a selected station for each sample size')
#ggsave(paste0(figdir,'stations_selected_by_multiobj_opt_m500.jpeg'))


#### Is there a correlation between near-surface lapse rates based on these samples and free-air lapse rates?
  source('Code/get_free_air_lapse_rates.R')
  free_air <- get_free_air_lapse_rates(sub_bbox,rtse=T,meta)
  
  # join free_air df with seasons to make it compatible with the station lapse rate calculations:
  seasons_numeric <- c(rep(c('Spring','Summer','Fall','Winter'),each=3),rep('Annual',12))
  seasons_numeric <- data.frame(cbind(seasons_numeric, c(3,4,5,6,7,8,9,10,11,12,1,2,1:12)))
  names(seasons_numeric) <- c('season','month')
  free_air <- join(free_air,seasons_numeric, by='month')
  # calculate the seasonal mean lapse rates as the average of the monthly lapse rates already calculated
  free_air <- free_air %>% group_by(year,lat,lon,season) %>% summarise(lapse_seasonal = mean(lapse, na.rm=T))
  # average seasonal lapse rates over the domain for comparison with stations:
  free_air <- free_air %>% group_by(year,season) %>% summarise(lapse=mean(lapse_seasonal, na.rm=T))
  
  
  # for each season and each sample size, find the correlation of free-air and near-surface lapse rates:
  facor <- data.frame('season'=rep(seas,each=(total_stations-1)), 'nn'= rep(2:total_stations,length(seas)), 'R'=NA)
  yrs <- c(1980:2016)
  for (ss in 1:length(seas)){
    fa <- free_air %>% filter(season==seas[ss], year %in% yrs)
    mins <- samp_wt_wide %>% filter(season==seas[ss]) %>% group_by(num_stations) %>% summarise(vmin =which.min(vsum))
    for (nn in 2:total_stations){
      station_ids <- t(combn(as.character(unique(meta_sc$station_id)),nn))[mins$vmin[nn-1],]
      meta_trim <- meta %>% filter(station_id %in% station_ids,season==seas[ss], year %in% yrs)
      plot(meta_trim$tmax,meta_trim$elev)
      ns <- meta_trim %>% group_by(year) %>% summarise(lapse=summary(lm(tmax~elev))$coefficients[2]*1000)
      facor$R[which(facor$season==seas[ss] & facor$nn==nn)] <- cor(fa$lapse,ns$lapse)^2
      plot(ns$year,ns$lapse, ylim=c(-7,1), type='l', col='blue')
      lines(fa$year, fa$lapse)
    }
  }
  ggplot(facor) + geom_point(aes(x=nn, y=R, col=season))



#################









#########
dat <- samp_wt_wide %>% filter(season=='Annual')
plot(dat$vsum, dat$lapse)

library(colorRamps)
ggplot(dat, aes(x=elev, y=m500, z=srad)) + geom_point(aes(col=lapse)) +
  scale_color_gradientn(name='lapse rate\n(C/km)',colors=matlab.like(50))


library(plotly)
colors <- matlab.like(100)
p <- plot_ly(dat, x = ~elev, y = ~m500, z = ~srad, color = ~lapse, colors = colors)
p







