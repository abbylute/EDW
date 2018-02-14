# meta should contain metadata for the stations to plot, including 'lat', 'lon', and 'elev'
# brks is the latitudinal or longitudinal bands to aggregate elevation over.  Should include the min and max of the range as well.
# title is a title for the whole plot.

map_stations_transects <- function(meta, brks, cols='black', maintitle){
  
  source('Documents/IGERT/EDW/Temperature_analysis/Code/get_elev_transects.R')
  library(ggplot2)
  library(cowplot) #plot_grid
  
  if (sum(c('lat','lon','elev') %in% names(meta))!=3){
    return(print('Meta argument must contain metadata for stations to plot.  This should inlude lat, lon, and elev fields with these names'))
  }
  
  meta <- cbind(meta,cols)
  
  midbrks <- (brks[1:(length(brks)-1)] + brks[2:length(brks)])/2
  meta$latbin <- factor(cut(meta$lat, brks), levels= sort(unique(cut(meta$lat, brks)),decreasing=T),labels=rev(midbrks))
  
  banddir <- 'Documents/IGERT/EDW/Temperature_analysis/'
  if (file.exists(paste0(banddir,'bandmax_',paste(brks,collapse=''),'.txt'))){
    bandmax <- read.table(paste0(banddir,'bandmax_',paste(brks,collapse=''),'.txt'))
    bandmean <- read.table(paste0(banddir,'bandmean_',paste(brks,collapse=''),'.txt'))
    bandmin <- read.table(paste0(banddir,'bandmin_',paste(brks,collapse=''),'.txt'))
  } else {
    bandmax <- get_elev_transects(brks,orient='horizontal',fun=max)
    bandmean <- get_elev_transects(brks,orient='horizontal',fun=mean)
    bandmin <- get_elev_transects(brks,orient='horizontal',fun=min)
  }
  
  colrng <- range(cols)
  p <- list()
  for (i in 1:length(midbrks)){
    # subset the latitude band:
    bx <- bandmax[which(bandmax$facet == midbrks[i]),]
    bm <- bandmean[which(bandmean$facet == midbrks[i]),]
    bn <- bandmin[which(bandmin$facet == midbrks[i]),]
    st <- meta[which(meta$latbin == midbrks[i]),]
    # plot:
   p[[i]] <- ggplot(bx, aes(x = lon, y = elev)) + geom_line(color='blue') +
             geom_line(inherit.aes=F, data=bm, aes(x=lon, y=elev), color='green') +
             geom_line(inherit.aes=F, data=bn, aes(x=lon, y=elev), color='brown') +
             geom_point(inherit.aes=F, data=st, aes(x=lon,y=elev,color=cols), size=.6) + 
             scale_color_gradientn(colors=rainbow(50)[15:50],'Trend C/decade',limits=c(colrng)) +
             labs(ylab=as.character(midbrks[i])) + theme_bw() + ylab(as.character(midbrks[i])) +
             theme(axis.title.x=element_blank()) +
             coord_cartesian(xlim = c(-125, -100))
  }

title <- ggdraw() + draw_label(maintitle, fontface='bold')
figs <- do.call(plot_grid,c(rev(p),ncol=1))
pp <- plot_grid(title,figs, ncol=1, rel_heights=c(.1,1))

pp
}

