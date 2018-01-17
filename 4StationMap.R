# Make map of Peyto Glacier region showing station locations.

# load map packages:
mapdir <-'/Volumes/Macintosh HD/Users/abbylute/Documents/IGERT/Coldwater/AbbysProject/Code/Map/'
dir <- "Documents/IGERT/Hydroclimatology/ResearchProject/TopoWx_Station_Data/Rockies/" 
source(paste0(mapdir,'FIX_get_stamenmap_RCode.R')) #ignore png warning
source(paste0(mapdir,'getbreaks.R'))
source(paste0(mapdir,'colorbar2sided.R'))
library(maptools) # fortify
library(rgdal) # for readOGR
library(ggplot2) #fortify
library(ggmap) #get_stamenmap
library(foreign) # to open dbfs
library(plyr) #join
library(scales) # for map scales
library(grid) # for 'unit'
library(gridExtra) # for tableGrob
library(gtable) # for arranging plot components
library(rgeos) 
library(raster)
library(colorRamps)

# load data:
  tmax <- read.table(paste0(dir,'Tmax_trend.txt'))
  tmin <- read.table(paste0(dir,'Tmin_trend.txt'))
  

  
  
# get a basemap for your bounding box
  win <- 4.5
  mbbox <- c(-126, 33,-103,50)
  gm <- get_stamenmap(bbox=mbbox,maptype="toner",zoom=6)  
  ggmap(gm)
  bb <- attr(gm,"bb") #extract the bounding box of the map

  
  latlines <- c(35,38,41,44,47,50)
  df2 <- data.frame(matrix(NA,ncol=5,nrow=6))
  names(df2) <- c('bin','lat','lon','latend','lonend')
  df2$bin <- c(1:6)
  df2$lat <- latlines
  df2$lon <- rep(-103,6)
  df2$latend <- latlines
  df2$lonend <- rep(-126,6)
  
gg <- ggmap(gm) +
  geom_segment(inherit.aes=FALSE, data=df2, aes(x=lon, y=lat, xend=lonend, yend=latend),size=.8,lty=2,col='grey') + 
  geom_point(inherit.aes=FALSE, data=tmax, aes(x=Longitude, y=Latitude,color=Elevation),size=1)+
  scale_color_gradientn(name="Station\nElevation (m)", colors=terrain.colors(20))+
  theme(legend.key=element_rect(fill="white")) +
  coord_map(ylim=mbbox[c(2,4)], xlim=mbbox[c(1,3)]) +
  labs(y='',x='')
gg





