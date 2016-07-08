#============================================================================================
# h1 : peak - left (main)
# h2 : rite - peak (main)
# delta : rite - left = h1 + h2 (main) 
# sub.delta : rite - left (sub)
get_patternFeature <- function(data, sub.data, speed, patternIndex, c_factor, use_Ap = FALSE, useConsistencyFlag = TRUE){
  
  timestamp <- data$timestamp
  value     <- data$value
  sub.value <- sub.data$value
  
  patternFeature        <- patternIndex[,c('str','end')]
  patternFeature$len    <- patternIndex$end - patternIndex$str + 1
  patternFeature$str.t  <- timestamp[patternIndex$str  ]
  patternFeature$end.t  <- timestamp[patternIndex$end+1]
  
  if('maxSlope' %in% names(patternIndex)) patternFeature$max.slope <- speed[ patternIndex$maxSlope ]
  if('minSlope' %in% names(patternIndex)) patternFeature$min.slope <- speed[ patternIndex$minSlope ]
  
  if('peak' %in% names(patternIndex)){
    peak.value.set           <- value[ patternIndex$peak ]
    patternFeature$h1        <- peak.value.set - value[ patternIndex$str ]
    patternFeature$h2        <- value[ patternIndex$end+1 ] - peak.value.set
  }
  patternFeature$sub.delta <- sub.value[patternIndex$end+1] - sub.value[patternIndex$str]
  patternFeature$delta     <- value[patternIndex$end+1]     - value[patternIndex$str]
  
  if( !missing(c_factor) && c_factor >= 1 ){
    
    patternFeature$org_h2        <-  patternFeature$h2
    patternFeature$org_sub.delta <-  patternFeature$sub.delta
    
    patternFeature <- subset( patternFeature, (str > c_factor*len) & (end < length(value)-c_factor*len) )
    if( nrow(patternFeature) == 0 ) return(NULL)
    
    consistency_flag <- rep(1, nrow(patternFeature))
    if( !use_Ap ){
      consistency_flag[ which(patternFeature$org_h2 > 0 | patternFeature$org_sub.delta > 0) ] <- 0
    }else{
      consistency_flag[ which(patternFeature$org_h2 > 0) ] <- 0
    }
      
    for(c_idx in 1:c_factor) {
      
      left.str <- patternFeature$str - patternFeature$len * c_idx
      left.end <- patternFeature$str - patternFeature$len * (c_idx-1) - 1
      rite.str <- patternFeature$end + patternFeature$len * (c_idx-1) + 2
      rite.end <- patternFeature$end + patternFeature$len * c_idx + 1
      
      tmp_sign_delta <- mapply( function(ls,le,rs,re) mean(value[rs:re]) - mean(value[ls:le]),
                             left.str, left.end, rite.str, rite.end )
      
      tmp_sign_sub.delta <- mapply( function(ls,le,rs,re) mean(sub.value[rs:re]) - mean(sub.value[ls:le]), 
                                    left.str, left.end, rite.str, rite.end)
      if( use_Ap ) tmp_sign_sub.delta <- 0 * tmp_sign_sub.delta
      
      consistency_flag[ which(tmp_sign_delta > 0 & tmp_sign_sub.delta > 0) ] <- 0
      
      patternFeature$h2        <- patternFeature$h2        + tmp_sign_delta
      patternFeature$delta     <- patternFeature$delta     + tmp_sign_delta
      patternFeature$sub.delta <- patternFeature$sub.delta + tmp_sign_sub.delta
      
    }
    
    if( useConsistencyFlag ) patternFeature <- patternFeature[which(consistency_flag == 1),]
    patternFeature$h2        <-  patternFeature$h2        / (c_factor+1)
    patternFeature$delta     <-  patternFeature$delta     / (c_factor+1)
    patternFeature$sub.delta <-  patternFeature$sub.delta / (c_factor+1)  
  }
  names(patternFeature) <- mapvalues( names(patternFeature), warn_missing=F, 
                                      c('str','end','len','str.t','end.t'), 
                                      c('start.idx','end.idx','pattern.size',
                                        'start.timestamp','end.timestamp') )
  
  return(patternFeature)
}
