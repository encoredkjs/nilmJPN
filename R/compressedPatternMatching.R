compressedPatternMatching <- function( compressedState, targetPattern ){
  
  pattern.length <- length(targetPattern) 
  
  if( targetPattern[1] %in% c(1,2) ){ # rising edge
    
    if( any( targetPattern %in% c(3,4) & shift(targetPattern,-1) %in% c(1,2) ) ){
      # power가 증가 -> 감소 -> 증가로 정의되면 마지막 증가는 패턴에서 제외
      pattern.length <- max(which(targetPattern %in% c(3,4) & shift(targetPattern,-1) %in% c(1,2) ))
    }
    
    matchedPattern.str <- grepl.pattern( compressedState$val, targetPattern )
    matchedPattern.end <- matchedPattern.str + pattern.length - 1
    maxSlope.pos <- matchedPattern.str + which( targetPattern == 1 & shift(targetPattern,-1) %in% c(2,3)) - 1
    minSlope.pos <- matchedPattern.str + which( targetPattern == 3 & shift(targetPattern,-1) %in% c(4,1)) - 1 
    peak.pos <- matchedPattern.str + which( targetPattern %in% c(1,2) & shift(targetPattern,-1) %in% c(3,4)) - 1
    
    result <- data.frame(str = compressedState$str[ matchedPattern.str ], 
                         end = compressedState$end[ matchedPattern.end ])
    
    if( length(maxSlope.pos) != 0 ) result <- cbind( result, data.frame(maxSlope = compressedState$end[ maxSlope.pos ]))
    if( length(minSlope.pos) != 0 ) result <- cbind( result, data.frame(minSlope = compressedState$end[ minSlope.pos ]))
    if( length(peak.pos) == 0 ) peak.pos <- matchedPattern.end
    
    result <- cbind( result, data.frame(peak = compressedState$end[ peak.pos ] + 1)) 
    return( result )
    
  }else if( targetPattern[1] %in% c(3,4) ){ # falling edge 
    
    matchedPattern.str <- grepl.pattern( compressedState$val, targetPattern )
    matchedPattern.end <- matchedPattern.str + pattern.length - 1
    minSlope.pos <- matchedPattern.str + which( targetPattern == 3 & shift(targetPattern,-1) %in% c(4,1,2)) - 1 
    peak.pos <- matchedPattern.str 
    
    return( data.frame(str = compressedState$str[ matchedPattern.str ], 
                       end = compressedState$end[ matchedPattern.end ],
                       minSlope = compressedState$end[ minSlope.pos ],
                       peak = compressedState$str[ peak.pos ]))
  }
}
