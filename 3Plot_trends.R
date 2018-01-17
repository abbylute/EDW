# Looking at EDW in Rockies.
# Plot TopoWx station trends
#--------------------------------------------------------
library(cowplot) #plot_grid
library(ggplot2)

dir <-   dir <- 'Documents/IGERT/Hydroclimatology/ResearchProject/TopoWx_Station_Data/Rockies/'

tmax <- read.table(paste0(dir,'Tmax_trend.txt'))
tmin <- read.table(paste0(dir,'Tmin_trend.txt'))

#----------------------------------------------------
# Bin stations by elevation and latitude:
Elevation_bin <- cut(pmax(tmax$Elevation,200.1), breaks = seq(200,3600,200), labels =as.character(seq(200,3400,200)))
plot(tmax$Elevation,Elevation_bin)
Lat_bin <- cut(tmax$Latitude, breaks = seq(35,50,3), labels = c('35_38','38_41','41_44','44_47','47_50'))
plot(tmax$Latitude,Lat_bin)
tmax <- cbind(tmax,Elevation_bin,Lat_bin)


Elevation_bin <- cut(pmax(tmin$Elevation,200.1), breaks = seq(200,3600,200), labels =as.character(seq(200,3400,200)))
plot(tmin$Elevation,Elevation_bin)
Lat_bin <- cut(tmin$Latitude, breaks = seq(35,50,3), labels = c('35_38','38_41','41_44','44_47','47_50'))
plot(tmin$Latitude,Lat_bin)
tmin <- cbind(tmin,Elevation_bin,Lat_bin)

#there was a small bit of correlation between elevation and latitude within 5 degree latitude bins, so I am boing to use smaller bins:
summary(lm(tmax$Elevation[which(tmax$Lat_bin=='35_40')]~tmax$Latitude[which(tmax$Lat_bin=='35_40')]))
summary(lm(tmax$Elevation[which(tmax$Lat_bin=='40_45')]~tmax$Latitude[which(tmax$Lat_bin=='40_45')]))
summary(lm(tmax$Elevation[which(tmax$Lat_bin=='45_50')]~tmax$Latitude[which(tmax$Lat_bin=='45_50')]))

#----------------------------------------------------
#### Try rearranging to create a 5x2 plot for each season showing tmax and tmin trends for every latitude band
  seasons <- c('Annual','Winter','Spring','Summer','Fall')
  vars <- c('tmax','tmin')
  latbins <- rev(as.character(unique(tmax$Lat_bin)))
  pval <- 0.05
  ylims <- c(pmin(min(tmax$Elevation),min(tmin$Elevation)), pmax(max(tmax$Elevation),max(tmax$Elevation)))
  
  season <- 1
  var <- 1
  pname <- which(names(tmax)==paste0(seasons[season],'_p'))
  sname <- paste0(seasons[season],'_slope')
  xlims <- c(-.2,.4) #c(-.33,.6)
  
  # This version just shows the distribution of points within an elevation band as a line
    for (i in 1:length(latbins)){
      assign(paste0('p',i),ggplot(tmax[which(tmax$Lat_bin==latbins[i] & tmax[,pname]<pval),], aes(x = get(sname), y = as.numeric(Elevation_bin),group=as.numeric(Elevation_bin))) + 
        geom_boxplot() +theme_bw() + #geom_jitter() +
          scale_y_continuous(name = "Elevation (m)" ,limits=c(1,17.5),breaks=seq(1,18,2),labels=as.character(seq(200,3600,400))) +
          scale_x_continuous(name= "Trend (C/decade)",limits=xlims))
      assign(paste0('p',i+5),ggplot(tmin[which(tmin$Lat_bin==latbins[i] & tmin[,pname]<pval),], aes(x = get(sname), y = as.numeric(Elevation_bin),group=as.numeric(Elevation_bin))) + 
        geom_boxplot() +theme_bw() + #geom_jitter() +
          scale_y_continuous(name = "Elevation (m)" ,limits=c(1,17.5),breaks=seq(1,18,2),labels=as.character(seq(200,3600,400))) +
          scale_x_continuous(name= "Trend (C/decade)",limits=xlims))
    }
    
    jpeg(filename=paste0(dir,'Annual_slopes.jpeg'),units='in',width=6,height=10,quality=100,res=500)
    plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,  ncol = 2)
    dev.off()






# -----------------------------------------------------------------------------------------
# This version shows the distribution of points within an elevation band as a boxplot and has a trend line with labels
    seasons <- c('Annual','Winter','Spring','Summer','Fall')
    latbins <- rev(as.character(unique(tmax$Lat_bin)))
    pval <- 0.05
    ylims <- c(pmin(min(tmax$Elevation),min(tmin$Elevation)), pmax(max(tmax$Elevation),max(tmax$Elevation)))
    
    season <- 5
    pname <- which(names(tmax)==paste0(seasons[season],'_p'))
    sname <- which(names(tmax)==paste0(seasons[season],'_slope'))
    xlims <- c(-.2,.4) #c(-.33,.6)
    latlabels <- c('50°N - 47°N','47°N - 44°N','44°N - 41°N','41°N - 38°N','38°N - 35°N')

for (i in 1:length(latbins)){
  fit <- summary(lm(tmax[which(tmax$Lat_bin==latbins[i] & tmax[,pname]<pval),sname] ~ tmax$Elevation[which(tmax$Lat_bin==latbins[i] & tmax[,pname]<pval)]))
  if (fit$coefficients[8]<.01){
    plab <- 'p < .01'
  }else {
    plab <- paste("p == ",round(fit$coefficients[8],3))
  }
  assign(paste0('p',i),ggplot(tmax[which(tmax$Lat_bin==latbins[i] & tmax[,pname]<pval),], aes(x = Elevation, y =get(names(tmax)[sname]),group=as.numeric(Elevation_bin))) + 
           geom_boxplot(outlier.size=.5)+ coord_flip() +#geom_jitter()+
           geom_smooth(inherit.aes =F, aes(x=Elevation,y=get(names(tmax)[sname])), method = "lm",fullrange=T,level=0.95) +
           scale_x_continuous(name = "Elevation (m)" ,limits=c(200,3600),breaks=seq(200,3600,400),labels=as.character(seq(200,3600,400))) +
           scale_y_continuous(name= "",limits=xlims)+ theme_bw() + theme(plot.margin = unit(c(0.1, 0.2, -.3, 0.2), "cm")) +
           annotate("text", x =3500, y = -.15, size = 2, parse = TRUE, label = paste("R^2 == ",round(fit$r.squared,3))) +
           annotate("text", x =3200, y = -.15, size = 2, parse = TRUE, label = plab))
  
  fit <- summary(lm(tmin[which(tmin$Lat_bin==latbins[i] & tmin[,pname]<pval),sname] ~ tmin$Elevation[which(tmin$Lat_bin==latbins[i] & tmin[,pname]<pval)]))
  if (fit$coefficients[8]<.01){
    plab <- 'p < .01'
  }else {
    plab <- paste("p == ",round(fit$coefficients[8],3))
  }
  assign(paste0('p',i+5),ggplot(tmin[which(tmin$Lat_bin==latbins[i] & tmin[,pname]<pval),], aes(x = Elevation, y =get(names(tmax)[sname]),group=as.numeric(Elevation_bin))) + 
           geom_boxplot(outlier.size=.5)+ coord_flip() +#geom_jitter()+
           geom_smooth(inherit.aes =F, aes(x=Elevation,y=get(names(tmax)[sname])), method = "lm",fullrange=T,level=0.95) +
           scale_x_continuous(name = "",limits=c(200,3600),breaks=seq(200,3600,400),labels=as.character(seq(200,3600,400))) +
           scale_y_continuous(name="" ,limits=xlims)+ theme_bw() +theme(plot.margin = unit(c(0.1, 0, -.3, -.2), "cm")) +
           annotate("text", x =3500, y = -.15, size = 2, parse = TRUE, label = paste("R^2 == ",round(fit$r.squared,3))) +
           annotate("text", x =3200, y = -.15, size = 2, parse = TRUE, label = plab))
}
    # additional pieces to add to the plot:
px1 <- ggplot() + geom_text(aes(x=6,y=0),label="Trend (°C/decade)",size=4)+
  xlim(c(0,10)) +  scale_y_continuous(expand=c(0,0))+ 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
  axis.text.y=element_blank(),axis.ticks=element_blank(),
  axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",
  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),plot.background=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 
px2 <- px1 + xlim(c(0,11))
pblank <- ggplot() + geom_text(aes(x=6,y=0),label="",size=4)+
  xlim(c(0,10)) +  scale_y_continuous(expand=c(0,0))+ 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 
for (i in 1:length(latbins)){
  assign(paste0('py',i),ggplot() + geom_text(aes(x=0,y=0),label=latlabels[i],size=4,angle=270)+
  scale_x_continuous(expand=c(0,0)) +  scale_y_continuous(expand=c(0,0))+ 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())+theme(plot.margin = unit(c(0, 0.1, 0, 0), "cm")))
}
ptmax <- ggplot() + geom_text(aes(x=6,y=0),label=paste0(seasons[season]," Tmax"),size=6)+
  scale_y_continuous(expand=c(0,0))+ xlim(c(0,10))+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 
ptmin <- ggplot() + geom_text(aes(x=5,y=0),label=paste0(seasons[season]," Tmin"),size=6)+
  scale_y_continuous(expand=c(0,0))+ xlim(c(0,10)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 

jpeg(filename=paste0(dir,seasons[season],'_slopes_boxplot_trendline.jpeg'),units='in',width=6,height=10,quality=100,res=500)
plot_grid(ptmax,ptmin,pblank,
          p1,p6,py1,
          p2,p7,py2,
          p3,p8,py3,
          p4,p9,py4,
          p5,p10,py5,
          px1,px2,pblank,  ncol = 3, rel_heights=c(.5,rep(3,5),.3),rel_widths=c(3,3,.2))
dev.off()





### _____________________________________________________________
pval <- 0.05
tmaxp <- c(length(which(tmax$Annual_p<pval)),
           length(which(tmax$Winter_p<pval)),
           length(which(tmax$Spring_p<pval)),
           length(which(tmax$Summer_p<pval)),
           length(which(tmax$Fall_p  <pval)))
tminp <- c(length(which(tmin$Annual_p<pval)),
           length(which(tmin$Winter_p<pval)),
           length(which(tmin$Spring_p<pval)),
           length(which(tmin$Summer_p<pval)),
           length(which(tmin$Fall_p  <pval)))
df <- rbind(tmaxp,tminp)/1681*100
jpeg(filename=paste0(dir,'Number_sig_trends_by_season.jpeg'),units='in',width=6,height=5,quality=100,res=500)
barplot(df,beside=T,col=c('red','blue'),names.arg=c('Annual','Winter','Spring','Summer','Fall'),
        ylab='% significant trends (p<0.05)',legend.text=c('Tmax','Tmin'),ylim=c(0,100))
dev.off()


rsq <- data.frame(matrix(NA,ncol=5,nrow=5*5*2))
names(rsq) <- c('T','season','lat','rsq','pval')
rsq$T <- rep(c('Tmax','Tmin'),each=25)
rsq$season <- rep(rep(seasons, each=5),2)
rsq$lat <- rep(latbins,10)

for (s in 1:length(seasons)){
  for (i in 1:length(latbins)){
    pname <- which(names(tmax)==paste0(seasons[s],'_p'))
    sname <- which(names(tmax)==paste0(seasons[s],'_slope'))
    
    fit <- summary(lm(tmax[which(tmax$Lat_bin==latbins[i] & tmax[,pname]<pval),sname] ~ tmax$Elevation[which(tmax$Lat_bin==latbins[i] & tmax[,pname]<pval)]))
    rsq$rsq[which(rsq$season==seasons[s]& rsq$lat==latbins[i] & rsq$T=='Tmax')] <- fit$r.squared
    rsq$pval[which(rsq$season==seasons[s]& rsq$lat==latbins[i] & rsq$T=='Tmax')] <- fit$coefficients[8]
    
    fit <- summary(lm(tmin[which(tmin$Lat_bin==latbins[i] & tmin[,pname]<pval),sname] ~ tmin$Elevation[which(tmin$Lat_bin==latbins[i] & tmin[,pname]<pval)]))
    rsq$rsq[which(rsq$season==seasons[s]& rsq$lat==latbins[i] & rsq$T=='Tmin')] <- fit$r.squared
    rsq$pval[which(rsq$season==seasons[s]& rsq$lat==latbins[i] & rsq$T=='Tmin')] <- fit$coefficients[8]
    
  }
}
    
rtmax <- ggplot()+geom_tile(data=rsq[1:25,],aes(x=season,y=lat,fill=rsq))+ title('Tmax')
rtmin <- ggplot()+geom_tile(data=rsq[26:50,],aes(x=season,y=lat,fill=rsq))+ title('Tmin')

ptmax <- ggplot()+geom_tile(data=rsq[1:25,],aes(x=season,y=lat,fill=pval))
ptmin <- ggplot()+geom_tile(data=rsq[26:50,],aes(x=season,y=lat,fill=pval))

jpeg(filename=paste0(dir,'Rsquared_season_lat_raster.jpeg'),units='in',width=9,height=5,quality=100,res=500)
plot_grid(rtmax,rtmin)
dev.off()
    
write.csv(rsq,paste0(dir,'Rsq_table.csv'))
