# Identify SNOTEL temperature shift dates:

#library(ncdf4)
library(tidyverse)

#ncfile <- 'Data/TopoWx/stn_obs_tmax.nc'

#nc <- nc_open(ncfile,write=F)
#lats <- ncvar_get(nc, 'latitude')
#lons <- ncvar_get(nc, 'longitude')

ha <- read.csv('Data/TopoWx/homog_adjust.csv')
snotels <- ha[(which(substr(ha$STN_ID,1,4) == 'NRCS')),]


# find the earliest date when SNOTEL temperatures had to be adjusted by more than 0.1C.  Stations not in this table are assumed to have good data.
sno_tab <- snotels %>% group_by(STN_ID) %>% filter(ADJ.C. > .1) %>% summarise(brktime = min(YEAR_MONTH_START))
ggplot(sno_tab) + geom_point(aes(x=STN_ID,y=brktime))

write.table(sno_tab, 'Data/TopoWx/SNOTEL_shift_date.txt')
