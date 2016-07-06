extractBoxShape <- function( power, rEdgeStr, rEdgeEnd, fEdgeStr, fEdgeEnd, 
                             buffer.size = 14, threshold = 10, debug.mode = F, search.stableRegion = F
                             , recursiveBoxRemoving = T ){
  
  if( fEdgeStr <= rEdgeEnd ){
    print('"fEdgeStr" should be greater than "rEdgeEnd"')
  }
  
  n <- length(power)
  buffer.size <- 14
  left.ind <- pmax( rEdgeStr - buffer.size:0, 1 )
  rite.ind <- pmin( fEdgeEnd + 0:buffer.size, n )

  boxReverse <- power
  if( search.stableRegion ){
    stationary.left <- find.stationary.state( boxReverse[left.ind],forward=F)[[1]]
    stationary.rite <- find.stationary.state( boxReverse[rite.ind],forward=T)[[1]]
    stationary.left <- max( rEdgeStr - stationary.left, 1)
    stationary.rite <- min( fEdgeEnd + stationary.rite, n)
  }else{
    stationary.left <- rEdgeStr
    stationary.rite <- fEdgeEnd
  }
  
  subinterval.ind <- (rEdgeEnd+1):fEdgeStr
  subinterval <- boxReverse[subinterval.ind]
    
  base.line <- min( boxReverse[c(stationary.left,stationary.rite+1)] )
  if( !all(base.line >= subinterval) ){
    
    boxReverse[(stationary.left+1):stationary.rite] <- base.line
    
    if(recursiveBoxRemoving){
      while( length(which(abs(diff(subinterval)) > threshold)) > 0 ){
        
        subbox <- pmax( subinterval - find.box.shape.reverse( subinterval ), 0 )
        boxReverse[subinterval.ind] <- boxReverse[subinterval.ind] + subbox  
        subinterval <- subinterval - subbox
        
        subinterval.lower <- envelopeDetector( data.frame(active_power = subinterval,reactive_power=0), 
                                               use.active.power=T )$LowerEnvelope
        
        if(!all((subinterval-subinterval.lower)==0)){
          boxReverse[subinterval.ind] <- boxReverse[subinterval.ind] + subinterval - subinterval.lower
          subinterval <- subinterval.lower
        }
        
        if( all(subbox==0) ) break
      }
      
      rEdgeLen <- rEdgeEnd - stationary.left + 1
      if( rEdgeLen == 1 ){
        boxReverse[ stationary.left+1 ] <- boxReverse[stationary.left]
      }else if( rEdgeLen > 1 ){
        boxReverse[ stationary.left:rEdgeEnd ] <- seq( boxReverse[stationary.left], boxReverse[rEdgeEnd],
                                                       length.out=rEdgeLen)
      }
      
      fEdgeLen <- stationary.rite - fEdgeStr + 1
      if( fEdgeLen == 1 ){
        boxReverse[ stationary.rite-1 ] <- boxReverse[stationary.rite]
      }else if( fEdgeLen > 1 ){
        boxReverse[ fEdgeStr:stationary.rite ] <- seq( boxReverse[fEdgeStr], boxReverse[stationary.rite],
                                                       length.out=fEdgeLen)
      }
    }
    boxReverse[ boxReverse > power ] <- power[ boxReverse > power ]    
  }
  
  boxReverse[ boxReverse > power ] <- power[ boxReverse > power ]
  boxShape <- power - boxReverse
  
  box.lists <- series.to.box.lists(1:length(boxShape), boxShape, 0 )[[1]]
  if( is.null(box.lists) ) return(power)
  
  box.lists <- subset( box.lists, duration > 5 )
  if( is.null(box.lists) | is.null(nrow(box.lists)) | nrow(box.lists) == 0 ) return( power )
  
  boxReverse.idx <- unique( unlist( mapply( function(str,end) str:end, c(1,box.lists$end+1), c(box.lists$str-1,n) )))
  boxReverse[boxReverse.idx] <- power[boxReverse.idx]
  boxShape <- power - boxReverse
  
  return( data.frame( 'boxShape' = boxShape, 'boxReverse' = boxReverse ) )
}
  