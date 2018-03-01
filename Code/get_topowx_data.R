
get_topowx_data <- function(varname, startdate, enddate,
                            bbox, exclude_network=character(), min_pdays, min_pyrs) {
  
  if (varname %in% c('tmax_raw','tmax_homog','tmax_infilled','tmin_raw','tmin_homog','tmin_infilled')){
    var <- substr(varname,1,4)
  }else{
    return(print('Error: varname argument must be one of: tmax_raw, tmax_homog, tmax_infilled, tmin_raw, tmin_homog, or tmin_infilled'))
  }

  library(ncdf4)
  library(mgcv)
  datadir <- 'Data/TopoWx/'
  #outdir <- 'Documents/IGERT/EDW/Temperature_analysis/'
  nc <- nc_open(paste0(datadir,'stn_obs_',var,'.nc'),write=FALSE)

  tt <- ncvar_get(nc,'time')  # insert dim names shown in print out to display dimension values
  #ncatt_get(nc,'time','units')
  cal <- seq(as.Date('1948-01-01'),as.Date('2016-12-31'),'days') #should change this to be more robust!!!!!!!!!
  if (!is.null(which(cal==startdate))){
    st <- which(cal == startdate)
  } else{
     return(print(paste0('startdate is not within timeframe, must be between ',cal[1],' and ',cal[length(cal)])))
  }
  if (!is.null(which(cal==enddate))){
    en <- which(cal == enddate)
  } else {
    return(print(paste0('enddate is not within timeframe, must be between ',cal[1],' and ',cal[length(cal)])))
  }

  lats <- ncvar_get(nc,'latitude')
  lons <- ncvar_get(nc,'longitude')
  elev <- ncvar_get(nc,'elevation')
  station_id <- as.character(ncvar_get(nc,'station_id'))
  network <- substr(station_id,1,4)
  snotel_meta <- read.csv('Data/SNOTEL/SNOTEL_metadata.csv',header=T,skip=988)
  xx <- which(substr(station_id,6,nchar(station_id)) %in% snotel_meta$Acton.Id)
  network[xx] <- "SNOTEL"
  
  # find stations within the bounding box (bbox)
  poly <- matrix(c(rep(bbox[1],2),rep(bbox[3],2),bbox[1],bbox[c(2,4,4,2)],bbox[2]), ncol=2)
  ins <- in.out(poly, matrix(c(lons,lats),ncol=2))
  #plot(poly,type='l')
  #points(lons[which(ins==TRUE)],lats[which(ins==TRUE)])


  temp <- ncvar_get(nc,start=c(1,st),count=c(-1,length(cal[st:en])),varname) # dims= station, day
  nc_close(nc)
  cal <- cal[st:en]
  
  # Trim dataset based on bounding box and excluded networks:
  if(length(exclude_network)==0){
    temp <- temp[which(ins==T),]
    elev <- elev[which(ins==T)]
    lat <- lats[which(ins==T)]
    lon <- lons[which(ins==T)]
    station_id <- station_id[which(ins==T)]
    network <- network[which(ins==T)]
  } else {
    temp <- temp[(which(ins==T & !network == exclude_network)),]
    elev <- elev[(which(ins==T & !network == exclude_network))]
    lat <- lats[(which(ins==T & !network == exclude_network))]
    lon <- lons[(which(ins==T & !network == exclude_network))]
    station_id <- station_id[(which(ins==T & !network == exclude_network))]
    network <- network[(which(ins==T & !network == exclude_network))]
  }

  # Trim dataset based on completeness of data:
  nstations <- nrow(temp)
  ny <- round(length(cal)/365)
  oct1 <- which(substr(cal,6,10)=='10-01')
  sep30 <- which(substr(cal,6,10)=='09-30')
  which_yrs_complete <- array(data=NA, dim = c(nstations,ny)) # stations, years, tmax type
  
  for (yy in 1:ny){ # cycle through years
    which_yrs_complete[,yy] <- apply(temp[,oct1[yy]:sep30[yy]], 1, function(x) (length(which(!is.na(x)))/length(x)>=min_pdays))
  }
  # identify the number of years with x% complete data at each stations since some date:
  nyrs_complete <- rowSums(which_yrs_complete)

  complete <- which((nyrs_complete/ny) > min_pyrs)
  temp <- temp[complete,]

  temp1 <- data.frame(cbind(station_id[complete],network[complete],lat[complete],lon[complete],
                            elev[complete]), stringsAsFactors=F)
  names(temp1) <- c('station_id','network','lat','lon','elev')

  out <- list(temp1,temp)
  print(paste0('This function returns a list with two items: the first (out[[1]]) is the metadata, 
              the second (out[[2]]) is the matrix of temperature values (rows=stations, cols=times)'))
out


} # end function
