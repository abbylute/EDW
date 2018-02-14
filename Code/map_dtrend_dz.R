
# trend is a dataframe containing temperature trend values for each station
# grid_res is the resolution of the grid of dtrend/dz evaluation points in degrees lat/lon (e.g. 1 means latitude spacing of 1 degree and longitude spacing of 1 degree)
# search_radius is in degrees lat/lon and is the radius around each grid point within which stations are incorporated for dtrend/dz calculations
# min_stations is the minimum number of stations within the search_radius of a grid point that is required to compute dtrend/dz
# plot_title is the title that will be put on the resulting plot
# bbox is the bounding box within which computations are completed

map_dtrend_dz <- function(trend, var, grid_res, search_radius, min_stations, figdir, figname, bbox=c(-125,30,-100,53)){

  library(ggplot2)
  library(cowplot) #plot_grid
  library(colorRamps)
  library(dplyr)
  library(ggmap) #get_stamenmap
  seasons <- c('Annual','Winter','Spring','Summer','Fall')
  gm <- get_stamenmap(bbox,maptype="toner",zoom=6)  
  
  
  for (season in 1:length(seasons)){
    dat <- select(trend, lat, lon, elev, paste0(seasons[season],'_slope'))
    grid_lats <- seq(bbox[2],bbox[4],grid_res)
    grid_lons <- seq(bbox[1],bbox[3],grid_res)
    griddf <- data.frame(lats=rep(grid_lats,length(grid_lons)), lons=rep(grid_lons,each=length(grid_lats)), dtdz=NA)

    # for each grid point (row of griddf), find the stations within 1 degree radius and compute dtdz
    for (ii in 1:nrow(griddf)){
      glat <- griddf$lat[ii]
      glon <- griddf$lon[ii]
      st <- dat[dat$lat>(glat-search_radius) & dat$lat<(glat+search_radius) & dat$lon>(glon-search_radius) & dat$lon<(glon+search_radius),]
      mod <- summary(lm(st[[paste0(seasons[season],'_slope')]] ~ st$elev))
      
      if (nrow(st)>=min_stations & diff(range(st$elev))>500 & mod$coefficients[8]<pval){
        griddf$dtdz[ii] <- mod$coefficients[2] *1000
      }
    }
  
    collim <- max(abs(range(griddf$dtdz,na.rm=T)))
    collim <- c(-collim,collim)
    gg <- ggmap(gm) + geom_tile(data=griddf, aes(x=lons, y=lats, fill=dtdz)) +
      scale_fill_gradientn('dTrend/dz\n(C/decade/km)',colors=matlab.like(20),na.value="transparent",limits=collim) +
      labs(x='',y='',title=paste0(seasons[season],' ',var,' EDW'))
    
    jpeg(filename=paste0(figdir,figname,'_',seasons[season],'_',var,'_dtdz_map.jpeg'),width=6,height=7,units="in",res=500,quality=100)
    print(gg)
    dev.off()
  print(seasons[season])
  }
}



# An earlier attempt with some cool use of new tidy skills:
#library(tidyr)
#gridmat <- as.matrix(arrange(spread(griddf,lons,dtdz),desc(lats)))[,-1]
#image(gridmat,  axes=FALSE)
#image(x=grid_lons,y=grid_lats,z=t(gridmat))
