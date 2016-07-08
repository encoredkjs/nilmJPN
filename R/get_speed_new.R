get_speed_new <- function(timestamp, value, allowInfSlope = F){
  
  if(class(timestamp)[1] != 'numeric') timestamp <- as.numeric(timestamp)
  time.diff  <- diff(timestamp)
  value.diff <- diff(value)
  slope      <- value.diff / time.diff
  
  if( !allowInfSlope ){
    infinite.slope <- which(!is.finite(slope))
    if( length(infinite.slope) > 0 ){
      finite.slope   <- setdiff( 1:length(slope), infinite.slope )
      if(             1 %in% infinite.slope ) slope[1:min(finite.slope)] <- slope[min(infinite.slope)]
      if( length(slope) %in% infinite.slope ) slope[max(finite.slope):length(slope)] <- slope[max(finite.slope)]
      infinite.slope <- setdiff( infinite.slope, c(1,length(slope)) )
#       idx <-   sapply( infinite.slope, function(x) c( max(finite.slope[ finite.slope < x ]), 
#                                                        min(finite.slope[ finite.slope > x ]) ))
#       slope[ infinite.slope ] <- ( slope[idx[1,]] + slope[idx[2,]] )/2
      split.infSlope <- split( infinite.slope, cumsum(c(-1,diff(infinite.slope)) != 1) )
      idx      <- sapply( split.infSlope, function(x) c( head(x,1)-1, tail(x,1)+1) )
      idx.left <- unlist(mapply( function(minVal,duplicatedItem) rep(minVal,length(duplicatedItem)), 
                                 idx[1,], split.infSlope ))
      idx.rite <- unlist(mapply( function(maxVal,duplicatedItem) rep(maxVal,length(duplicatedItem)), 
                                 idx[2,], split.infSlope ))
      slope[ infinite.slope ] <- ( slope[idx.left] + slope[idx.rite] )/2
    }
  }
  if( length(which(!is.finite(slope))) > 0 ) stop("Algorithm is not correct : get_speed_new")
  return( slope )
}
