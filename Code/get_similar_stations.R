# This function identifies similar stations of a specified sample size based 
# on similarity of specified parameters (e.g. topographic position, solar radiation, etc)

# meta: a dataframe containing a column of identifiers (e.g. station id) and a column for each of
# the specified 'vars' that will be used to calculate similarity

# vars: a string or vector of strings of variable names to use in computing similarity (e.g. c('srad','tpim500'))

# opt: a vector the length of vars indicating whether the range of each var should be minimized or maximized. Example: opt=c('min','max') for vars=c('tpi','elev')

# num_stations: the number of stations to output

# scale: if scale=T (default), then normalize the vars

get_similar_stations = function(meta, vars, opt, num_stations, scale=T){
  if (!(hasArg(opt) & hasArg(meta) & hasArg(vars) & hasArg(num_stations))){
    print('Error: argument(s) missing.  Arguments must include meta, vars, opt, and num_stations.  scale arg is optional.')
  }
  if (!all(c('station_id',vars) %in% names(meta))){
    print('Error: dataframe called meta must contain column names station_id and the names of the specified variables to use for calculating similarity (vars)')
  }
  
  if (num_stations>nrow(meta)){
    print('Error: the number of stations to return is greater than the number of stations in the supplied dataset (meta)')
  }
  
  if (!length(which(opt %in% c('max','min')))==length(vars)){
    print('Error: the opt argument must be the same length as the vars argument and must contain only max or min strings')
  }
  
  opt_lookup <- data.frame('max'=-1, 'min'=1)
  
  cc <- which(names(meta) %in% vars)
  
  if (scale==T){
    meta[,cc] <- scale(meta[,cc])
  }
  
  samp <- combn(meta$station_id,num_stations) # compute all possible combinations of this number of stations
  
  samp_vars <- matrix(nrow=ncol(samp), ncol=length(vars)) # set up an table with the ranges of each variable for each sample
  
  # compute the ranges of each variable for each sample:
  for (ii in 1:ncol(samp)) {
    tt <- meta[which(meta$station_id %in% samp[,ii]),]
    for (vv in 1:length(vars)){
      samp_vars[ii,vv] <- diff(range(tt[,cc[vv]])) # because I'm using the range function, there's no need to use the abs value, it does this (indirectly) automatically
    }
  }
  
  # invert ranges that need to be maximized:
  opts <- opt_lookup[opt]
  samp_vars <- apply(X=samp_vars,oo=opts, FUN = function(X,oo){out <- X^oo; return(out)}, MARGIN=1)

  stations <- samp[,which.min(sqrt(colSums(samp_vars^2)))] # this sample minimizes the variability of the vars
  #stations <- samp[,which.min(rowSums(samp_vars))] the above method seems like it balances the multiple vars better than this method
  return(stations)
}