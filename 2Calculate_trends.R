# Looking at EDW in Rockies.
# Calculate trends and significance
#--------------------------------------------------------

# Set Up:
dir <-   dir <- 'Documents/IGERT/Hydroclimatology/ResearchProject/TopoWx_Station_Data/Rockies/'

tmax <- as.matrix(read.table(paste0(dir,'Tmax.txt')))
tmaxmeta <- read.table(paste0(dir,'Tmax_meta.txt'))
tmin <- as.matrix(read.table(paste0(dir,'Tmin.txt')))
tminmeta <- read.table(paste0(dir,'Tmin_meta.txt'))

cal <- seq(as.Date('1948-01-01'),as.Date('2015-12-31'),'days')

#-----------------------------------------------------------------------
# calculate annual and seasonal Tmax and Tmin values:

styr <- which(substr(cal,6,10)=='01-01'); styr <- c(styr,ncol(tmax))
aTmax <- matrix(NA,nrow=nrow(tmax),ncol=(length(styr)-1))
for (ss in 1:nrow(tmax)){
  for(yy in 1:(length(styr)-1)){
    aTmax[ss,yy] <- mean(tmax[ss,styr[yy]:(styr[yy+1]-1)])
  }
}

styr <- which(substr(cal,6,10)=='01-01'); styr <- c(styr,ncol(tmin))
aTmin <- matrix(NA,nrow=nrow(tmin),ncol=(length(styr)-1))
for (ss in 1:nrow(tmin)){
  for(yy in 1:(length(styr)-1)){
    aTmin[ss,yy] <- mean(tmin[ss,styr[yy]:(styr[yy+1]-1)])
  }
}

Tmax_summer <- matrix(NA,nrow=nrow(tmax),ncol=(length(styr)-1)) # June, July, August
Tmax_fall   <- matrix(NA,nrow=nrow(tmax),ncol=(length(styr)-1)) # Sept, Oct, Nov
Tmax_winter <- matrix(NA,nrow=nrow(tmax),ncol=(length(styr)-1)) # Dec, Jan, Feb
Tmax_spring <- matrix(NA,nrow=nrow(tmax),ncol=(length(styr)-1)) # Mar, Apr, May
for (ss in 1:nrow(tmax)){
  for(yy in 1:(length(styr)-1)){
    Tmax_summer[ss,yy] <- mean(tmax[ss,(styr[yy]+which(substr(cal[styr[yy]:(styr[yy+1]-1)],6,7) %in% c('06','07','08'))-1)])
    Tmax_fall[ss,yy]   <- mean(tmax[ss,(styr[yy]+which(substr(cal[styr[yy]:(styr[yy+1]-1)],6,7) %in% c('09','10','11'))-1)])
    Tmax_winter[ss,yy] <- mean(tmax[ss,(styr[yy]+which(substr(cal[styr[yy]:(styr[yy+1]-1)],6,7) %in% c('12','01','02'))-1)])
    Tmax_spring[ss,yy] <- mean(tmax[ss,(styr[yy]+which(substr(cal[styr[yy]:(styr[yy+1]-1)],6,7) %in% c('03','04','05'))-1)])
  }
}

Tmin_summer <- matrix(NA,nrow=nrow(tmin),ncol=(length(styr)-1)) # June, July, August
Tmin_fall <- matrix(NA,nrow=nrow(tmin),ncol=(length(styr)-1)) # Sept, Oct, Nov
Tmin_winter <- matrix(NA,nrow=nrow(tmin),ncol=(length(styr)-1)) # Dec, Jan, Feb
Tmin_spring <- matrix(NA,nrow=nrow(tmin),ncol=(length(styr)-1)) # Mar, Apr, May
for (ss in 1:nrow(tmin)){
  for(yy in 1:(length(styr)-1)){
    Tmin_summer[ss,yy] <- mean(tmin[ss,(styr[yy]+which(substr(cal[styr[yy]:(styr[yy+1]-1)],6,7) %in% c('06','07','08'))-1)])
    Tmin_fall[ss,yy] <- mean(tmin[ss,(styr[yy]+which(substr(cal[styr[yy]:(styr[yy+1]-1)],6,7) %in% c('09','10','11'))-1)])
    Tmin_winter[ss,yy] <- mean(tmin[ss,(styr[yy]+which(substr(cal[styr[yy]:(styr[yy+1]-1)],6,7) %in% c('12','01','02'))-1)])
    Tmin_spring[ss,yy] <- mean(tmin[ss,(styr[yy]+which(substr(cal[styr[yy]:(styr[yy+1]-1)],6,7) %in% c('03','04','05'))-1)])
  }
}



#-----------------------------------------------------------------------
# annual and seasonal trends and signficance (using sen's slope and MK test)
library(trend)

Tmax_trend <- data.frame(matrix(NA,nrow=nrow(tmax),ncol=10))
Tmin_trend <- data.frame(matrix(NA,nrow=nrow(tmin),ncol=10))
Tmax_trend <- cbind(tmaxmeta,Tmax_trend)
Tmin_trend <- cbind(tminmeta,Tmin_trend)
names(Tmax_trend) <- c('Elevation','Latitude','Longitude','Annual_slope','Annual_p','Winter_slope','Winter_p','Spring_slope','Spring_p','Summer_slope','Summer_p','Fall_slope','Fall_p')
names(Tmin_trend) <- names(Tmax_trend)


Tmax_trend$Annual_slope <- apply(aTmax,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
Tmax_trend$Winter_slope <- apply(Tmax_winter,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
Tmax_trend$Spring_slope <- apply(Tmax_spring,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
Tmax_trend$Summer_slope <- apply(Tmax_summer,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
Tmax_trend$Fall_slope   <- apply(Tmax_fall,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)

Tmax_trend$Annual_p <- apply(aTmax,1,function(x) mk.test(x)$p.value)
Tmax_trend$Winter_p <- apply(Tmax_winter,1,function(x) mk.test(x)$p.value)
Tmax_trend$Spring_p <- apply(Tmax_spring,1,function(x) mk.test(x)$p.value)
Tmax_trend$Summer_p <- apply(Tmax_summer,1,function(x) mk.test(x)$p.value)
Tmax_trend$Fall_p   <- apply(Tmax_fall,1,function(x) mk.test(x)$p.value)


Tmin_trend$Annual_slope <- apply(aTmin,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
Tmin_trend$Winter_slope <- apply(Tmin_winter,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
Tmin_trend$Spring_slope <- apply(Tmin_spring,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
Tmin_trend$Summer_slope <- apply(Tmin_summer,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
Tmin_trend$Fall_slope <- apply(Tmin_fall,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)

Tmin_trend$Annual_p <- apply(aTmin,1,function(x) mk.test(x)$p.value)
Tmin_trend$Winter_p <- apply(Tmin_winter,1,function(x) mk.test(x)$p.value)
Tmin_trend$Spring_p <- apply(Tmin_spring,1,function(x) mk.test(x)$p.value)
Tmin_trend$Summer_p <- apply(Tmin_summer,1,function(x) mk.test(x)$p.value)
Tmin_trend$Fall_p <- apply(Tmin_fall,1,function(x) mk.test(x)$p.value)

write.table(Tmax_trend,paste0(dir,'Tmax_trend.txt'))
write.table(Tmin_trend,paste0(dir,'Tmin_trend.txt'))








# old way:
sslopep <- apply(aTmax,1,function(x) sens.slope(x,0.95)$p.value)
sslopep_summer <- apply(Tmax_summer,1,function(x) sens.slope(x,0.95)$p.value)
sslopep_fall <- apply(Tmax_fall,1,function(x) sens.slope(x,0.95)$p.value)
sslopep_winter <- apply(Tmax_winter,1,function(x) sens.slope(x,0.95)$p.value)
sslopep_spring <- apply(Tmax_spring,1,function(x) sens.slope(x,0.95)$p.value)

# Look at Tmax trends across elevations:   
# hard to see any relationship between annual Tmax trend and elevation for whole US Rockies
plot(sslope,elevin)
cor.test(sslope,elevin,method="pearson")
plot(sslope_summer,elevin)
cor.test(sslope_summer,elevin,method="pearson")
plot(sslope_fall,elevin)
cor.test(sslope_fall,elevin,method="pearson")
plot(sslope_winter,elevin)
cor.test(sslope_winter,elevin,method="pearson")
plot(sslope_spring,elevin)
cor.test(sslope_spring,elevin,method="pearson")

# might be better to break this into a few geographic pieces though since the elevation values are conflating
# N-S Rockies heights (higher in S) and the amount of elevation contrast across a given slope
latin <- lats[which(ins==TRUE)]
lonsin <- lons[which(ins==TRUE)]
inN <- in.out(polyN, matrix(c(lonsin,latin),ncol=2))
inS <- in.out(polyS,matrix(c(lonsin,latin),ncol=2))

plot(sslope[inN],elevin[inN])
cor.test(sslope[inN],elevin[inN],method="pearson")

plot(sslope[inS],elevin[inS])
cor.test(sslope[inS],elevin[inS],method="pearson")

plot(sslope_summer[inN],elevin[inN])
cor.test(sslope_summer[inN],elevin[inN],method="pearson")
plot(sslope_fall[inN],elevin[inN])
cor.test(sslope_fall[inN],elevin[inN],method="pearson")
plot(sslope_winter[inN],elevin[inN])
cor.test(sslope_winter[inN],elevin[inN],method="pearson")
plot(sslope_spring[inN],elevin[inN])
cor.test(sslope_spring[inN],elevin[inN],method="pearson")
plot(sslope_summer[inS],elevin[inS])
cor.test(sslope_summer[inS],elevin[inS],method="pearson")
plot(sslope_fall[inS],elevin[inS])
cor.test(sslope_fall[inS],elevin[inS],method="pearson")
plot(sslope_winter[inS],elevin[inS])
cor.test(sslope_winter[inS],elevin[inS],method="pearson")
plot(sslope_spring[inS],elevin[inS])
cor.test(sslope_spring[inS],elevin[inS],method="pearson")

#elevink <- elevin/1000

alim <- c(-.25,.6)
filenm <- paste0('Documents/IGERT/Hydroclimatology/ResearchProject/Tmax_TopoWx_NSRockies_scatter.jpeg')
jpeg(filename=filenm,width=180,height=260,units="mm",res=3000,quality=100)
nf <- layout(mat = matrix(c(1:24),6,4,byrow=F), widths = c(.5,2,2,2), height = c(.5,rep(2,5)), TRUE)
layout.show(nf)
par(mar=c(1.8,2,0.1,0.2),tcl=-0.2, mgp=c(.8, .1, 0))

plot(0,0,xaxt='n',yaxt='n',bty='n',cex=0,ann=F)

plot(0,0,xaxt='n',yaxt='n',bty='n',cex=0,ann=F)
text(0,-.1,'Annual',cex=1.5,srt=90)

plot(0,0,xaxt='n',yaxt='n',bty='n',cex=0,ann=F)
text(0,-.1,'Summer',cex=1.5,srt=90)

plot(0,0,xaxt='n',yaxt='n',bty='n',cex=0,ann=F)
text(0,-.1,'Fall',cex=1.5,srt=90)

plot(0,0,xaxt='n',yaxt='n',bty='n',cex=0,ann=F)
text(0,-.1,'Winter',cex=1.5,srt=90)

plot(0,0,xaxt='n',yaxt='n',bty='n',cex=0,ann=F)
text(0,-.1,'Spring',cex=1.5,srt=90)

plot(0,0,xaxt='n',yaxt='n',bty='n',cex=0,ann=F)
text(0,-.1,'Rocky Mountains',cex=1.5)

plot(sslope,elevin,xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin~(sslope)))
plot(sslope_summer,elevin,xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin~(sslope_summer)))
plot(sslope_fall,elevin,xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin~(sslope_fall)))
plot(sslope_winter,elevin,xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin~(sslope_winter)))
plot(sslope_spring,elevin,xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin~(sslope_spring)))

plot(0,0,xaxt='n',yaxt='n',bty='n',cex=0,ann=F)
text(0,-.1,'Northern Rockies',cex=1.5)

plot(sslope[inN],elevin[inN],xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin[inN]~(sslope[inN])))
plot(sslope_summer[inN],elevin[inN],xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin[inN]~(sslope_summer[inN])))
plot(sslope_fall[inN],elevin[inN],xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin[inN]~(sslope_fall[inN])))
plot(sslope_winter[inN],elevin[inN],xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin[inN]~(sslope_winter[inN])))
plot(sslope_spring[inN],elevin[inN],xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin[inN]~(sslope_spring[inN])))

plot(0,0,xaxt='n',yaxt='n',bty='n',cex=0,ann=F)
text(0,-.1,'Southern Rockies',cex=1.5)

plot(sslope[inS],elevin[inS],xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin[inS]~(sslope[inS])))
plot(sslope_summer[inS],elevin[inS],xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin[inS]~(sslope[inS])))
plot(sslope_fall[inS],elevin[inS],xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin[inS]~(sslope[inS])))
plot(sslope_winter[inS],elevin[inS],xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin[inS]~(sslope[inS])))
plot(sslope_spring[inS],elevin[inS],xlab='trend (C/decade)',ylab='elevation (m)',cex=.3,xlim=alim)
lines(c(0,0),c(0,4000),col='grey',lty=3)
abline(lm(elevin[inS]~(sslope[inS])))

dev.off()





###########  NOT COMPLETED BELOW THIS LINE ##############

#write.table(tmax,'Documents/IGERT/Hydroclimatology/ResearchProject/TopoWx_Data/TopoWx_Tmax_infilled_all.txt')
nc_close(fff)

write.table(cal,'Documents/IGERT/Hydroclimatology/ResearchProject/TopoWx_Data/TopoWx_calendar.txt')

# TMIN
fff <- paste0(dir,'stn_obs_tmin.nc')
nc <- nc_open(fff,write=FALSE)
print(nc)

tt1 <- ncvar_get(nc,'time')  # insert dim names shown in print out to display dimension values
identical(tt,tt1)
ntt <- dim(tt1) # 1948-2015
head(tt1)
ncatt_get(nc,'time','units')
#cal <- seq(as.Date('1948-01-01'),as.Date('2015-12-31'),'days')


lats <- ncvar_get(nc,'latitude')
lons <- ncvar_get(nc,'longitude')
elev <- ncvar_get(nc,'elevation')


dname <- "tmin_infilled"  # [stn_id,time]  # degrees C
tmax <- ncvar_get(nc,start=c(1,1),count=c(-1,-1),dname) # dims= station, day
dim(tmax)


filval <- 9.96920996838687e+36
tmax[which(tmax> (filval-.1))]
min(tmax)
max(tmax)
# appears there are not any missing values.

write.table(tmax,'Documents/IGERT/Hydroclimatology/ResearchProject/TopoWx_Tmax_infilled_all.txt')
nc_close(fff)



