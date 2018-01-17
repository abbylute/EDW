# Look for EDW in TopoWx station data for Rocky Mountain Region

# define region
  library(mgcv)
  x <- c(-119.5, -113, -104, -105, -109)
  y <- c(49, 49, 40, 35, 35)
  poly <-matrix(c(x,y),ncol=2)
  plot(x,y,type='l')
  
  #northern (ID, MT, WY)
  x <- c(-119.5, -113, -107.5, -117)
  y <- c(49, 49, 43, 43)
  polyN <-matrix(c(x,y),ncol=2)
  
  #southern (mostly CO, a bit of UT and NM)
  x <- c( -107.5, -117, -104, -105, -109)
  y <- c(43, 43, 40, 35, 35)
  polyS <-matrix(c(x,y),ncol=2)
  
# Step 1: extract temperature data for all sites (infilled version) for both Tmax and Tmin
  library(ncdf4)
  #dir <- '/Volumes/Seagate Expansion Drive/SNOTEL_Analogs/Data/TopoWx_Station_Data/'
  dir <- 'Documents/IGERT/Hydroclimatology/ResearchProject/TopoWx_Station_Data/'

# TMAX 
  fff <- paste0(dir,'stn_obs_tmax.nc')
  nc <- nc_open(fff,write=FALSE)
  print(nc)
  
  tt <- ncvar_get(nc,'time')  # insert dim names shown in print out to display dimension values
  ncatt_get(nc,'time','units')
  cal <- seq(as.Date('1948-01-01'),as.Date('2015-12-31'),'days')
  
  lats <- ncvar_get(nc,'latitude')
  lons <- ncvar_get(nc,'longitude')
  elev <- ncvar_get(nc,'elevation')

  ins <- in.out(poly, matrix(c(lons,lats),ncol=2))
  plot(poly,type='l')
  points(lons[which(ins==TRUE)],lats[which(ins==TRUE)])

  dname <- "tmax_infilled"  # [stn_id,time]  # degrees C
  tmax <- ncvar_get(nc,start=c(1,1),count=c(-1,-1),dname) # dims= station, day
  dim(tmax)
  tmaxin <- tmax[which(ins==TRUE),]
  dim(tmaxin)

  elevin <- elev[which(ins==TRUE)]
  latin <- lats[which(ins==TRUE)]
  lonin <- lons[which(ins==TRUE)]
  write.table(tmaxin, paste0(dir,'Rockies/Tmax.txt'))
  write.table(cbind(elevin,latin,lonin), paste0(dir,'Rockies/Tmax_meta.txt'))
  write.table(cal, paste0(dir,'Rockies/Tmax_cal.txt'))
  
  nc_close(nc)
  
################################
  # Repeat for Tmin 
  fff <- paste0(dir,'stn_obs_tmin.nc')
  nc <- nc_open(fff,write=FALSE)
  print(nc)
  
  tt <- ncvar_get(nc,'time')  # insert dim names shown in print out to display dimension values
  ncatt_get(nc,'time','units')
  cal <- seq(as.Date('1948-01-01'),as.Date('2015-12-31'),'days')
  
  lats <- ncvar_get(nc,'latitude')
  lons <- ncvar_get(nc,'longitude')
  elev <- ncvar_get(nc,'elevation')
  
  ins <- in.out(poly, matrix(c(lons,lats),ncol=2))
  plot(poly,type='l')
  points(lons[which(ins==TRUE)],lats[which(ins==TRUE)])
  
  dname <- "tmin_infilled"  # [stn_id,time]  # degrees C
  tmin <- ncvar_get(nc,start=c(1,1),count=c(-1,-1),dname) # dims= station, day
  dim(tmin)
  tminin <- tmin[which(ins==TRUE),]
  dim(tminin)
  
  elevin <- elev[which(ins==TRUE)]
  latin <- lats[which(ins==TRUE)]
  lonin <- lons[which(ins==TRUE)]
  write.table(tminin, paste0(dir,'Rockies/Tmin.txt'))
  write.table(cbind(elevin,latin,lonin), paste0(dir,'Rockies/Tmin_meta.txt'))
  write.table(cal, paste0(dir,'Rockies/Tmin_cal.txt'))
  
  nc_close(nc)
  
  
  
  
  
  
  
 