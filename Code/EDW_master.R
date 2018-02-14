# master script:

# 0. Set these variables:
#.........................................................
datatype <- 'homog' # 'raw','homog',or 'infilled'
varname1 <- paste0('tmax_',datatype)
varname2 <- paste0('tmin_',datatype)
startdate <- '1980-01-01' # these should correspond to calendar years to work with the trend analysis
enddate <- '2015-12-31'
exclude <-  character()#'SNOTEL' # character()#default is
bbox <- c(-125,30,-100,53)
min_pdays <- 0.85  # minimum percent of days of data within each year
min_pyrs <- 0.90 # minimum percent of years of data at each station for the period defined above
savefigs <- T
figname <- paste0('ALL_',datatype,'_',as.character(min_pdays),'pd_',as.character(min_pyrs),'py_st',substr(startdate,1,4))
figdir <- paste0('Figures/ALL_',datatype,'_',substr(startdate,1,4),'/')
brks <- c(30,34,38,42,46,53) # latitude breaks to create bins


codedir <- 'Code/'
seasons <- c('Annual','Winter','Spring','Summer','Fall')

# 1. Get TopoWx data
#.......................................................................
source(paste0(codedir,'get_topowx_data.R'))
tmax <- get_topowx_data(varname=varname1, startdate=startdate, enddate=enddate,
                        bbox=bbox, exclude_network=exclude, min_pdays=min_pdays, min_pyrs=min_pyrs)
tmax_meta <- tmax[[1]]
tmax_meta[,3:5] <- as.numeric(unlist(tmax_meta[,3:5]))
tmax <- tmax[[2]]

tmin <- get_topowx_data(varname=varname2, startdate=startdate, enddate=enddate,
                        bbox=bbox, exclude_network=exclude, min_pdays=min_pdays, min_pyrs=min_pyrs)
tmin_meta <- tmin[[1]]
tmin_meta[,3:5] <- as.numeric(unlist(tmin_meta[,3:5]))
tmin <- tmin[[2]]

# 2. Map the stations you selected
#.......................................................................
source('Code/FIX_get_stamenmap_RCode.R') #ignore png warning
library(ggmap) #get_stamenmap

gm <- get_stamenmap(bbox,maptype="toner",zoom=6)  

tmax_map <- ggmap(gm) +
  geom_point(data=tmax_meta, aes(x=lon, y=lat,color=elev),size=1)+
  scale_color_gradientn(name="Elevation", colors=rainbow(20))+
  theme(legend.key=element_rect(fill="white")) +
  coord_map(ylim=bbox[c(2,4)], xlim=bbox[c(1,3)]) +
  labs(y='', x='', title=paste0('TopoWx ',varname1,' stations (n=',nrow(tmax_meta),')\nwith at least ',min_pyrs*100,'% of years \nwith at least ',min_pdays*100,'% of data since ',substr(startdate,1,4)) )

tmin_map <- ggmap(gm) +
  geom_point(data=tmin_meta, aes(x=lon, y=lat,color=elev),size=1)+
  scale_color_gradientn(name="Elevation", colors=rainbow(20))+
  theme(legend.key=element_rect(fill="white")) +
  coord_map(ylim=bbox[c(2,4)], xlim=bbox[c(1,3)]) +
  labs(y='', x='', title=paste0('TopoWx ',varname2,' stations (n=',nrow(tmin_meta),')\nwith at least ',min_pyrs*100,'% of years \nwith at least ',min_pdays*100,'% of data since ',substr(startdate,1,4)) )

if (savefigs==T){
  jpeg(filename=paste0(figdir,figname,'_tmax_map.jpeg'),width=6,height=7,units="in",res=500,quality=100)
  print(tmax_map)
  dev.off()
  jpeg(filename=paste0(figdir,figname,'_tmin_map.jpeg'),width=6,height=7,units="in",res=500,quality=100)
  print(tmin_map)
  dev.off()
} else {
  print(tmax_map)
  print(tmin_map)
}

# 3. Plot the elevation transects with stations on them
#.......................................................................
source(paste0(codedir,'map_stations_transects.R'))

tmax_bandplot <- map_stations_transects(tmax_meta, brks, 
        title=paste0('TopoWx ',varname1,' stations (n=',nrow(tmax_meta),')\nwith at least ',min_pyrs*100,'% of years \nwith at least ',min_pdays*100,'% of data since ',substr(startdate,1,4)))
if (savefigs==T){
  jpeg(filename=paste0(figdir,figname,'_elev_lat_bandplot.jpeg'),width=6,height=7,units="in",res=500,quality=100)
  tmax_bandplot
  dev.off()
} else {
  tmax_bandplot
}

tmin_bandplot <- map_stations_transects(tmin_meta, brks, 
        title=paste0('TopoWx ',varname2,' stations (n=',nrow(tmin_meta),')\nwith at least ',min_pyrs*100,'% of years \nwith at least ',min_pdays*100,'% of data since ',substr(startdate,1,4)))
if (savefigs==T){
  jpeg(filename=paste0(figdir,figname,'_elev_lat_bandplot.jpeg'),width=6,height=7,units="in",res=500,quality=100)
  tmin_bandplot
  dev.off()
} else {
  tmin_bandplot
}





# 4. Calculate Trends
#.......................................................................
source(paste0(codedir,'get_trends.R'))
tmax_trends <- get_trends(data=tmax, meta=tmax_meta, startdate=startdate, enddate=enddate)
# need to adjust this function to accomodate missing values starting in line 47.
tmin_trends <- get_trends(data=tmin, meta=tmin_meta, startdate=startdate, enddate=enddate)



# 5. Plot Trends on Elevation Bands
#.......................................................................
source(paste0(codedir,'plot_EDW_trend_transects.R'))

for (season in 1:length(seasons)){
  cc <- which(names(tmax_trends)==paste0(seasons[season],'_slope'))
  tmax_bandplot <- plot_EDW_trend_transects(tmax_meta, brks, cols=tmax_trends[,cc],
                       maintitle=paste0('TopoWx ',seasons[season],' ',varname1,' stations (n=',nrow(tmax_meta),')\nwith at least ',min_pyrs*100,'% of years \nwith at least ',min_pdays*100,'% of data since ',substr(startdate,1,4)))
  tmin_bandplot <- plot_EDW_trend_transects(tmin_meta, brks, cols=tmin_trends[,cc],
                       maintitle=paste0('TopoWx ',seasons[season],' ',varname2,' stations (n=',nrow(tmin_meta),')\nwith at least ',min_pyrs*100,'% of years \nwith at least ',min_pdays*100,'% of data since ',substr(startdate,1,4)))
  if (savefigs==T){
    jpeg(filename=paste0(figdir,figname,'_',seasons[season],'_tmax_elev_lat_bandplot.jpeg'),width=6,height=8,units="in",res=500,quality=100)
    print(tmax_bandplot)
    dev.off()
    jpeg(filename=paste0(figdir,figname,'_',seasons[season],'_tmin_elev_lat_bandplot.jpeg'),width=6,height=8,units="in",res=500,quality=100)
    print(tmin_bandplot)
    dev.off()
  } else {
    print(tmax_bandplot)
    print(tmin_bandplot)
  }
}


# 5. Lattice Plot of Trends vs Elevation for each latitude bin and tmax/tmin
#..............................................................................
source(paste0(codedir,'make_EDW_trends_lattice_plot.R'))
    make_EDW_trends_lattice_plot(tmax_trends, tmin_trends, brks, varnames=c(varname1,varname2), pval=0.10, figdir)


# 5. Map Station Trends
#..............................................................................
    library(colorRamps)
    collims <- max(abs(range(tmax_trends[,c(6,8,10,12,14)], tmin_trends[,c(6,8,10,12,14)])))
    collims <- c(-collims, collims)
  for (season in 1:length(seasons)){
    col <- which(names(tmax_trends)==paste0(seasons[season],'_slope'))
    siz <- which(names(tmax_trends)==paste0(seasons[season],'_p'))

    tmax_map <- ggmap(gm) +
      geom_point(data=tmax_trends, aes_string(x='lon', y='lat',color=colnames(tmax_trends)[col],size=(((colnames(tmax_trends)[siz])))))+
      scale_color_gradientn(name="Trend C/decade", colors=matlab.like(40), limits=c(collims))+
      scale_size_continuous(range = c(1.6,.5), breaks=0.10, guide=F) +
      theme(legend.key=element_rect(fill="white")) +
      coord_map(ylim=bbox[c(2,4)], xlim=bbox[c(1,3)]) +
      labs(y='', x='', title=paste0('TopoWx ',seasons[season],' ',varname1,' stations (n=',nrow(tmax_meta),')\nwith at least ',min_pyrs*100,'% of years \nwith at least ',min_pdays*100,'% of data since ',substr(startdate,1,4)) )
    
    tmin_map <- ggmap(gm) +
      geom_point(data=tmin_trends, aes_string(x='lon', y='lat',color=colnames(tmin_trends)[col],size=(((colnames(tmin_trends)[siz])))))+
      scale_color_gradientn(name="Trend C/decade", colors=matlab.like(40), limits=c(collims))+
      scale_size_continuous(range = c(1.6,.5), breaks=0.10, guide=F) +
      theme(legend.key=element_rect(fill="white")) +
      coord_map(ylim=bbox[c(2,4)], xlim=bbox[c(1,3)]) +
      labs(y='', x='', title=paste0('TopoWx ',seasons[season],' ',varname2,' stations (n=',nrow(tmin_meta),')\nwith at least ',min_pyrs*100,'% of years \nwith at least ',min_pdays*100,'% of data since ',substr(startdate,1,4)) )
    
    if (savefigs==T){
      jpeg(filename=paste0(figdir,figname,'_',seasons[season],'_tmax_trend_map.jpeg'),width=6,height=7,units="in",res=500,quality=100)
      print(tmax_map)
      dev.off()
      jpeg(filename=paste0(figdir,figname,'_',seasons[season],'_tmin_trend_map.jpeg'),width=6,height=7,units="in",res=500,quality=100)
      print(tmin_map)
      dev.off()
    } else {
      print(tmax_map)
      print(tmin_map)
    }
  }
