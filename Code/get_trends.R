
# Calculate trends and significance
#--------------------------------------------------------
# df should be a dataframe (rows=stations, columns=timesteps)
# meta should be corresponding metadata (rows=stations)
# startdate and enddate define the time frame over which to compute trends (e.g. '1980-10-01' and '2015-09-30')

# at the moment this function returns annual and seasonal trends and significance.
# in the future it might be useful to modify this to compute trends over a specified (subannual) time frame


get_trends <- function(data, meta, startdate, enddate){

#-----------------------------------------------------------------------
# calculate annual and seasonal Tmax and Tmin values:
  cal <- seq(as.Date(startdate),as.Date(enddate),'days')

  styr <- which(substr(cal,6,10)=='01-01'); styr <- c(styr,(ncol(tmax)+1))  # This should be run over calendar years, not
  # water years.  If run over WY, selection of 'fall' dates below gets split across two years.
  years <- substr(cal[styr],1,4)
  annualdf <- matrix(NA,nrow=nrow(data),ncol=(length(styr)-1))
  for (ss in 1:nrow(data)){
    for(yy in 1:(length(styr)-1)){
      annualdf[ss,yy] <- mean(data[ss,styr[yy]:(styr[yy+1]-1)],na.rm=T)
    }
  }
  
  stfall <- which(substr(cal,6,10)=='09-01')
  stwinter <- which(substr(cal,6,10)=='12-01')
  stspring <- which(substr(cal,6,10)=='03-01')
  stsummer <- which(substr(cal,6,10)=='06-01')
  summerdf <- matrix(NA,nrow=nrow(data),ncol=(length(styr)-1)) # June, July, August
  falldf   <- matrix(NA,nrow=nrow(data),ncol=(length(styr)-1)) # Sept, Oct, Nov
  winterdf <- matrix(NA,nrow=nrow(data),ncol=(length(styr)-1)) # Dec, Jan, Feb
  springdf <- matrix(NA,nrow=nrow(data),ncol=(length(styr)-1)) # Mar, Apr, May
  for (ss in 1:nrow(data)){
    for(yy in 1:(length(styr)-1)){
      #summerdf[ss,yy] <- mean(data[ss,(styr[yy]+which(substr(cal[styr[yy]:(styr[yy+1]-1)],6,7) %in% c('06','07','08'))-1)],na.rm=T)
      #falldf[ss,yy]   <- mean(data[ss,(styr[yy]+which(substr(cal[styr[yy]:(styr[yy+1]-1)],6,7) %in% c('09','10','11'))-1)],na.rm=T)
      #winterdf[ss,yy] <- mean(data[ss,(styr[yy]+which(substr(cal[styr[yy]:(styr[yy+1]-1)],6,7) %in% c('12','01','02'))-1)],na.rm=T)
      #springdf[ss,yy] <- mean(data[ss,(styr[yy]+which(substr(cal[styr[yy]:(styr[yy+1]-1)],6,7) %in% c('03','04','05'))-1)],na.rm=T)
      
      springdf[ss,yy] <- mean(data[ss, stspring[yy]:(stsummer[yy]-1)], na.rm=T)
      summerdf[ss,yy] <- mean(data[ss, stsummer[yy]:(stfall[yy]-1)], na.rm=T)
      falldf[ss,yy] <- mean(data[ss, stfall[yy]:(stwinter[yy]-1)], na.rm=T)
      if (yy==1){
        winterdf[ss,yy] <- mean(data[ss, 1:(stspring[yy]-1)], na.rm=T)
      } else {
        winterdf[ss,yy] <- mean(data[ss, stwinter[yy-1]:(stspring[yy]-1)], na.rm=T)
      }
    }
  }

  
#-----------------------------------------------------------------------
# annual and seasonal trends and signficance (using sen's slope and MK test)
 " library(trend)
  trends <- data.frame(matrix(NA,nrow=nrow(data),ncol=10))
  trends <- cbind(meta,trends)
  names(trends) <- c(names(meta),'Annual_slope','Annual_p','Winter_slope','Winter_p','Spring_slope','Spring_p','Summer_slope','Summer_p','Fall_slope','Fall_p')
  
  trends$Annual_slope <- apply(annualdf,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
  trends$Winter_slope <- apply(winterdf,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
  trends$Spring_slope <- apply(springdf,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
  trends$Summer_slope <- apply(summerdf,1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
  trends$Fall_slope   <- apply(falldf,  1,function(x) sens.slope(x,0.95)$estimates)*10 #(C/Decade)
  
  trends$Annual_p <- apply(annualdf,1,function(x) mk.test(x)$p.value)
  trends$Winter_p <- apply(winterdf,1,function(x) mk.test(x)$p.value)
  trends$Spring_p <- apply(springdf,1,function(x) mk.test(x)$p.value)
  trends$Summer_p <- apply(summerdf,1,function(x) mk.test(x)$p.value)
  trends$Fall_p   <- apply(falldf,  1,function(x) mk.test(x)$p.value)"  # old version of trend test that didn't allow missing values


  # new version of trend test that does allow missing values
  library(rkt)
  trends <- data.frame(matrix(NA,nrow=nrow(data),ncol=10))
  trends <- cbind(meta,trends)
  names(trends) <- c(names(meta),'Annual_slope','Annual_p','Winter_slope','Winter_p','Spring_slope','Spring_p','Summer_slope','Summer_p','Fall_slope','Fall_p')
  
  years <- as.numeric(unique(substr(cal,1,4)))
  #rrr <- rkt(years, annualdf[1,])
  #sens.slope(annualdf[1,],0.95)$estimates
  #mk.test(annualdf[1,])$p.value
  
  trends$Annual_slope <- apply(annualdf,1,function(x) rkt(years,x)$B)*10 #(C/Decade)
  trends$Winter_slope <- apply(winterdf,1,function(x) rkt(years,x)$B)*10 #(C/Decade)
  trends$Spring_slope <- apply(springdf,1,function(x) rkt(years,x)$B)*10 #(C/Decade)
  trends$Summer_slope <- apply(summerdf,1,function(x) rkt(years,x)$B)*10 #(C/Decade)
  trends$Fall_slope   <- apply(falldf,  1,function(x) rkt(years,x)$B)*10 #(C/Decade)
  
  trends$Annual_p <- apply(annualdf,1,function(x) rkt(years,x)$sl)
  trends$Winter_p <- apply(winterdf,1,function(x) rkt(years,x)$sl)
  trends$Spring_p <- apply(springdf,1,function(x) rkt(years,x)$sl)
  trends$Summer_p <- apply(summerdf,1,function(x) rkt(years,x)$sl)
  trends$Fall_p   <- apply(falldf,  1,function(x) rkt(years,x)$sl)
  
#write.table(trends,paste0(dir,'trends.txt'))
  print('For temperatures, trends are in units of C/decade.')
  
  trends

} # end function