
make_EDW_trends_lattice_plot <- function(tmax_trends, tmin_trends, latbrks, varnames, pval, savedir) {
  library(cowplot) #plot_grid
  library(ggplot2)
  t1 <- tmax_trends
  t2 <- tmin_trends
  if(!length(varnames)==2){
    return(print('Length of varnames argument should be two but it is not. Varnames should be the name of the tmax variable followed by the name of the tmin variable.'))
  }

  if (all(!c('lat','lon','elev') %in% names(t1)) || all(!c('lat','lon','elev') %in% names(t2))){
    return(print('Tmax_trends and Tmin_trends dataframes must contain latitude, longitude, and elevation data with column names lat, lon, and elev'))
  }
  
  nl <- length(latbrks)-1
  latbins <- paste0(latbrks[1:(nl)],'_',latbrks[2:(nl+1)])

#----------------------------------------------------
  # Bin stations by elevation and latitude:
  Elevation_bin <- cut(pmax(t1$elev,200.1), breaks = seq(200,3600,200), labels =as.character(seq(200,3400,200)))
  Lat_bin <- cut(t1$lat, breaks = latbrks, labels = latbins)
  t1 <- cbind(t1,Elevation_bin,Lat_bin)
  
  Elevation_bin <- cut(pmax(t2$elev,200.1), breaks = seq(200,3600,200), labels =as.character(seq(200,3400,200)))
  Lat_bin <- cut(t2$lat, breaks = latbrks, labels = latbins)
  t2 <- cbind(t2,Elevation_bin,Lat_bin)
  

# -----------------------------------------------------------------------------------------
# This version shows the distribution of points within an elevation band as a boxplot and has a trend line with labels
  seasons <- c('Annual','Winter','Spring','Summer','Fall')
  #latbins <- rev(as.character(unique(tmax$Lat_bin)))
  ylims <- c(pmin(min(t1$elev),min(t2$elev)), pmax(max(t1$elev),max(t2$elev)))
  xlims <- c(-.2,.4) #c(-.33,.6)
  latlabels <- paste0(rev(substr(latbins,4,5)),'°N - ',rev(substr(latbins,1,2)),'°N')
  
# extra plot pieces:
  theme_blankplots <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                      axis.text.y=element_blank(),axis.ticks=element_blank(),
                      axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",
                      panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),plot.background=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 
  xlab1 <- ggplot() + geom_text(aes(x=6,y=0),label="Trend (°C/decade)",size=4)+
        xlim(c(0,10)) +  scale_y_continuous(expand=c(0,0))+ theme_blankplots
  xlab2 <- xlab1 + xlim(c(0,11)) #for tmin
  pblank <- ggplot() + geom_text(aes(x=6,y=0),label="",size=4)+
      xlim(c(0,10)) +  scale_y_continuous(expand=c(0,0))+ theme_blankplots
  latlabs <- list()
  for (i in seq_along(latbins)){
    #assign(paste0('py',i),
    latlabs[[i]] <- ggplot() + geom_text(aes(x=0,y=0),label=latlabels[i],size=4,angle=270)+
      scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0,0))+ theme_blankplots
  }


plots <- list()

for (season in 1:5){
  print(seasons[season])
  title1 <- ggplot() + geom_text(aes(x=6,y=0),label=paste0(seasons[season]," ",varnames[1]),size=6)+
    scale_y_continuous(expand=c(0,0))+ xlim(c(0,10))+ theme_blankplots
  title2 <- ggplot() + geom_text(aes(x=6,y=0),label=paste0(seasons[season]," ",varnames[2]),size=6)+
    scale_y_continuous(expand=c(0,0))+ xlim(c(0,10))+ theme_blankplots
  
  for (i in rev(seq_along(latbins))){ # this needs to be in reverse (with rev) so that plots go from N to S.
    # for tmax:
    pname <- which(names(t1)==paste0(seasons[season],'_p'))
    sname <- which(names(t1)==paste0(seasons[season],'_slope'))
    
    if (length(which(t1$Lat_bin==latbins[i] & t1[,pname]<pval))<3){
      plots[[(2*i)-1]] <- pblank
    } else {
      fit <- summary(lm(t1[which(t1$Lat_bin==latbins[i] & t1[,pname]<pval),sname] ~ t1$elev[which(t1$Lat_bin==latbins[i] & t1[,pname]<pval)]))
      if (fit$coefficients[8]<.01){
        plab <- 'p < .01'
      }else {
        plab <- paste("p == ",round(fit$coefficients[8],3))
      }
      #assign(paste0('max',i),
            plots[[(2*i)-1]]<- ggplot(t1[which(t1$Lat_bin==latbins[i] & t1[,pname]<pval),], aes(x = elev, y =get(names(t1)[sname]),group=as.numeric(Elevation_bin))) + 
               geom_boxplot(outlier.size=.5)+ coord_flip() +#geom_jitter()+
               geom_smooth(inherit.aes =F, aes(x=elev,y=get(names(t1)[sname])), method = "lm",fullrange=T,level=0.95) +
               scale_x_continuous(name = "Elevation (m)" ,limits=c(200,3600),breaks=seq(200,3600,400),labels=as.character(seq(200,3600,400))) +
               scale_y_continuous(name= "",limits=xlims)+ theme_bw() + theme(plot.margin = unit(c(0.1, 0.2, -.3, 0.2), "cm")) +
               annotate("text", x =3500, y = -.15, size = 2, parse = TRUE, label = paste("R^2 == ",round(fit$r.squared,3))) +
               annotate("text", x =3200, y = -.15, size = 2, parse = TRUE, label = plab)
    }     
    
    # for tmin:
    pname <- which(names(t2)==paste0(seasons[season],'_p'))
    sname <- which(names(t2)==paste0(seasons[season],'_slope'))

    if (length(which(t2$Lat_bin==latbins[i] & t2[,pname]<pval))<3){
      plots[[2*i]] <- pblank
    } else {
      fit <- summary(lm(t2[which(t2$Lat_bin==latbins[i] & t2[,pname]<pval),sname] ~ t2$elev[which(t2$Lat_bin==latbins[i] & t2[,pname]<pval)]))
      if (fit$coefficients[8]<.01){
        plab <- 'p < .01'
      }else {
        plab <- paste("p == ",round(fit$coefficients[8],3))
      }
      #assign(paste0('max',i),
      plots[[2*i]]<- ggplot(t2[which(t2$Lat_bin==latbins[i] & t2[,pname]<pval),], aes(x = elev, y =get(names(t2)[sname]),group=as.numeric(Elevation_bin))) + 
        geom_boxplot(outlier.size=.5)+ coord_flip() +#geom_jitter()+
        geom_smooth(inherit.aes =F, aes(x=elev,y=get(names(t2)[sname])), method = "lm",fullrange=T,level=0.95) +
        scale_x_continuous(name = "Elevation (m)" ,limits=c(200,3600),breaks=seq(200,3600,400),labels=as.character(seq(200,3600,400))) +
        scale_y_continuous(name= "",limits=xlims)+ theme_bw() + theme(plot.margin = unit(c(0.1, 0.2, -.3, 0.2), "cm")) +
        annotate("text", x =3500, y = -.15, size = 2, parse = TRUE, label = paste("R^2 == ",round(fit$r.squared,3))) +
        annotate("text", x =3200, y = -.15, size = 2, parse = TRUE, label = plab)
    }
  } # end latbins
  
  # save the plot:
  titles <- plot_grid(title1,title2,pblank,ncol=3,rel_widths=c(1,1,.3))
  ylabs <- do.call(plot_grid,c(latlabs,ncol=1))
  body <- do.call(plot_grid,c(plots,ncol=2))
  xlabs <- plot_grid(xlab1, xlab2, pblank, ncol=3, rel_widths=c(1,1,.3))
  masterplot <-plot_grid(titles, plot_grid(body,ylabs,ncol=2,rel_widths=c(1.3,.1)), xlabs, nrow=3, rel_heights=c(.1,1.5,.1))
  
  jpeg(filename=paste0(savedir,seasons[season],'_slopes_boxplot_trendline.jpeg'),units='in',width=6,height=10,quality=100,res=500)
  print(masterplot)
  dev.off()
  
} # end season
#plot(t2[which(t2$Lat_bin==latbins[i]),pname])
} # end function