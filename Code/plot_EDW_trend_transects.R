# meta should contain metadata for the stations to plot, including 'lat', 'lon', and 'elev'
# brks is the latitudinal or longitudinal bands to aggregate elevation over.  Should include the min and max of the range as well.
# title is a title for the whole plot.

plot_EDW_trend_transects <- function(meta, brks, cols='black', maintitle){
  
  source('Documents/IGERT/EDW/Temperature_analysis/Code/get_elev_transects.R')
  library(ggplot2)
  library(cowplot) #plot_grid
  library(colorRamps)
  
  if (sum(c('lat','lon','elev') %in% names(meta))!=3){
    return(print('Meta argument must contain metadata for stations to plot.  This should inlude lat, lon, and elev fields with these names'))
  }
  
  meta <- cbind(meta,cols)
  
  midbrks <- (brks[1:(length(brks)-1)] + brks[2:length(brks)])/2
  meta$latbin <- factor(cut(meta$lat, brks), levels= sort(unique(cut(meta$lat, brks)),decreasing=T),labels=rev(midbrks))
  colrng <- c(-max(abs(range(cols))), max(abs(range(cols))))
  p <- list()
  for (i in 1:length(midbrks)){

    st <- meta[which(meta$latbin == midbrks[i]),]
    # plot:
   p[[i]] <- ggplot(st, aes(x = lon, y = elev,color=cols))  +
             geom_point(size=1) + 
             scale_color_gradientn(colors=matlab.like(10),'Trend C/decade',limits=c(colrng)) +
             labs(ylab=as.character(midbrks[i])) + theme_bw() + ylab(as.character(midbrks[i])) +
             theme(axis.title.x=element_blank()) +
             coord_cartesian(xlim = c(-125, -100))
  }

title <- ggdraw() + draw_label(maintitle, fontface='bold')
figs <- do.call(plot_grid,c(rev(p),ncol=1))
pp <- plot_grid(title,figs, ncol=1, rel_heights=c(.1,1))

pp
}

